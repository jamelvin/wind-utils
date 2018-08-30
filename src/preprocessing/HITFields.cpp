//  Copyright 2016 National Renewhite Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applichite law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#include "HITFields.h"
#include "core/YamlUtils.h"
#include "core/KokkosWrappers.h"
#include "core/PerfUtils.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

#include <fstream>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, HITFields, "init_hit_fields");

HITFields::HITFields(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void HITFields::load(const YAML::Node& node)
{
    // Setup mean velocity field
    auto mvel = node["mean_velocity"].as<std::vector<double>>();
    if (mvel.size() !=3)
        throw std::runtime_error("Invalid mean velocity field provided");
    mean_vel_ = mvel;

    // Process part info
    auto fluid_partnames = node["fluid_parts"].as<std::vector<std::string>>();
    fluid_parts_.resize(fluid_partnames.size());

    auto& meta = mesh_.meta();
    for(size_t i=0; i < fluid_partnames.size(); i++) {
        auto* part = meta.get_part(fluid_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    }

    // Get the HIT filename
    ti_ = node["desired_turb_intensity"].as<double>();
    // Get the mesh dimensions
    urms_ = node["urms_database"].as<double>();
    // Get the requested initial query time for the database
    t0_ = node["initialTime_database"].as<double>();
    // Get the requested number of x points desired in the 
    // fringe region
    nx_ = node["num_x_nodes_in_fringe"].as<int>();
    // Get the requested number of points to sample in each
    // direction in the y-z plane 
    ny_ = nz_ = node["num_sample_points"].as<double>();
    // Get the requested number of time slices to store in the exodus db
    numSteps_ = node["number_of_time_slices"].as<int>();
    // Get the requested dt between slices to store in the exodus db
    dtPhysical_ = node["physical_time_between_slices"].as<double>();

    // FIXME: Should be able to read the mesh and get these...
    // Get the lengths of the physical domain
    physLen_ = node["physical_domain_lengths"].as<std::vector<double>>();
    // Get the number of points in the input mesh in x and y
    numPts_ = node["number_of_input_mesh_points"].as<std::vector<int>>();
    // Get the minimum x coord in the physical domain
    minx_ = node["minimum_physical_x_coord"].as<double>();
}

void HITFields::initialize()
{
    const std::string timerName = "HITields::initialize";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    VectorFieldType& velocity = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");
    for (auto part: fluid_parts_) {
        stk::mesh::put_field(velocity, *part);
    }
    mesh_.add_output_field("velocity");

    // I have moved this here so I can output multiple "timesteps"
    // FIXME: There is almost certainly a better way to do this...
    fh_ = mesh_.stkio_.create_output_mesh("exoTesting.exo", stk::io::WRITE_RESULTS);
    for (auto fname: mesh_.output_fields_ ) {
        stk::mesh::FieldBase* fld = stk::mesh::get_field_by_name(fname, meta);
        if (fld != NULL) {
            mesh_.stkio_.add_field(fh_, *fld, fname);
        }
    }

    // Calculate additional quantities based on input
    uSweep_ = urms_ / ti_;

    // Ly_ and Lz_ are the size of the y and z directions in the JHTDB
    // that will be sampled.  If your physical domain is a box in y and z
    // these would both be equal to 2*pi as you would directly map onto
    // your domain.
    //  
    // TODO: Figure out what is best here when physical domain is not equal
    // in y and z ... This is Tony's recommendation to take largest length 
    // scale and only partially sample in the shorter dir
//    Ly_ = Lz_ = 2.0*PI;
    Ly_ = 2.0*PI * physLen_[1] / std::max(physLen_[1], physLen_[2]);
    Lz_ = 2.0*PI * physLen_[2] / std::max(physLen_[1], physLen_[2]);

//    tf_ = t0_ + Ly_ * nx_ / uSweep_ / ny_;
//    Lx_ = (tf_ - t0_) * uSweep_;
    filtWidth_ = Ly_ / ny_;
}

void HITFields::run()
{
    double x[nx_];
    double tol = 0.001;
    double dtLoc, time;

    // These are the number of gridpoints in each y-z plane
    // I am using this to allocate correct array sizes
    int ny = numPts_[1];
    int nz = numPts_[2];    
    double dx = physLen_[0]/(numPts_[0] - 1);
    double xmin = minx_;

    // FIXME: figure out a way to get dx from the mesh file...
    for (int i = 0; i < nx_; i++)
        x[i] = xmin + i*dx; 

    //dt = (nx_ == 1) ? 0.0 : (tf_ - t0_) / (nx_ - 1);

    // This is the time increment in the JHTDB associated with sampling
    // at a new location, i.e. for initial conditions
    dtLoc = 2.0*PI * dx / uSweep_ / std::max(physLen_[1], physLen_[2]);

    // Timing
    const std::string timerName = "HITields::run";
    auto timeMon = get_stopwatch(timerName);

    // Get meta data setup
    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();
    const int nDim = meta.spatial_dimension();

    // Get the velocity and coordinates fields from the mesh
    VectorFieldType* velocity = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");
    VectorFieldType* coordinates = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    // Select the mesh regions associated with fluid
    stk::mesh::Selector sel = stk::mesh::selectUnion(fluid_parts_);
    auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);

    // initialize vectors for JHTDB
    float pts[ny*nz][3] = {0};
    float vels[ny*nz][3] = {0};

    // set up the soap call to communicate with JHTDB
    init_JHTDB();

    // Loop over multiple "timesteps" to generate additional slices for the
    // inflow file
    
    for (int ts = 0; ts < numSteps_; ts++) {
      // This is the dt experienced by the database associated with adding a
      // slice at a physical time dtPhysical.  Note that this is negative
      // since it is like rewinding time back to where the slice would have
      // been to get to the boundary point at this physical time.
      double dtBdry = -dtPhysical_ * ts
                            * (2.0 * PI / std::max(physLen_[1], physLen_[2]))
                            * mean_vel_[0] / uSweep_;

      // Loop through each y-z plane, match the points in the physical domain
      // with corresponding pts in the database, query the velocities and
      // add them back into the mesh file
      //
      // Each y-z plane is queried from a particular x location and time in
      // the database.  As we move to the next plane, the x location and time
      // in the database are swept by the uSweep velocity, thus we need to 
      // query each separately
      for (int i = 0; i < nx_; i++) { 
  
        // We will loop through every grid point, and will use j as a counter
        // to determine the current position in the pts and vels array to add
        // the next eligible point for that y-z plane
        int j = 0;
  
        // The time for this y-z slice from the JHTDB
        // Note that dtBdry includes the multiplication by the loop variable
        // in its definition above.
        time = t0_ + dtLoc*i + dtBdry;
 
        // The x position for this y-z slice from the JHTDB
        double xcoord = uSweep_ * (dtLoc*i + dtBdry - t0_);
  
        // Loop through the mesh points
        for (size_t ib=0; ib < bkts.size(); ib++) {
          auto& b = *bkts[ib];
          for (size_t in=0; in < b.size(); in++) {
            auto node = b[in];
  
            // Retrieve the vel and coords at this node
            double* coords = stk::mesh::field_data(*coordinates, node);
            double* vel = stk::mesh::field_data(*velocity, node);
  
            // Only select points in this particular y-z plane of the mesh
            if (fabs(coords[0] - x[i]) < tol) {
              // Mapping of physical coordinates to JHTDB coordinates -- Tony's rec
              pts[j][0] = xcoord;
              pts[j][1] = coords[1] / std::max(physLen_[1], physLen_[2]) * 2.0 * PI;
              pts[j][2] = coords[2] / std::max(physLen_[1], physLen_[2]) * 2.0 * PI;
             
              // While we are here, lets initialize the velocities to 0
              vel[0] = 0.0;
              vel[1] = 0.0;
              vel[2] = 0.0;
   
              // Increment our array position
              j = j + 1;
            }
          }
        }

        std::cout << "Sampling JHTDB at x = " << xcoord << " and t = " << time << std::endl;
        if (j != ny*nz)
           std::cout << "WARNING: The number of grid points found for the y-z plane, "
                     << j << ", doesn't match ny*nz, " << ny*nz << ".  This seems like "
                     << "an error." << std::endl;
        
        // Query the JHTDB
        getBoxFilter (authtoken_.c_str(), dataset_.c_str(), field_.c_str(), 
                           time, filtWidth_, ny*nz, pts, vels);
  
        // Reset j for the next time through the arrays
        j = 0;
  
        // Loop through the mesh points
        for (size_t ib=0; ib < bkts.size(); ib++) {
          auto& b = *bkts[ib];
          for (size_t in=0; in < b.size(); in++) {
            auto node = b[in];
  
            // Retrieve the vel and coords at this node
            double* vel = stk::mesh::field_data(*velocity, node);
            double* coords = stk::mesh::field_data(*coordinates, node);
  
            // Only select points in this particular y-z plane of the mesh
            if (fabs(coords[0] - x[i]) < tol) {
              // The vels array is properly aligned with the order we will 
              // hit each node, so just directly copy the velocity in.  We
              // need to add the convective velocity here too as this will
              // be used in the real runs to overwrite the velocity in the
              // fringe region of the domain.
              vel[0] = vels[j][0] * mean_vel_[0]/uSweep_ + mean_vel_[0];
              vel[1] = vels[j][1] * mean_vel_[0]/uSweep_ + mean_vel_[1];
              vel[2] = vels[j][2] * mean_vel_[0]/uSweep_ + mean_vel_[2];
  
              // Increment our array position
              j = j + 1;
            }
          }
        }
      }
  
      // adding this to write multiple "timesteps"
      mesh_.stkio_.begin_output_step(fh_, dtPhysical_ * ts);
      mesh_.stkio_.write_defined_output_fields(fh_);
      mesh_.stkio_.end_output_step(fh_);
    }  // End timestep loop

    // Finalize the soap process for communication with JHTDB
    close_JHTDB();
}

void HITFields::init_JHTDB()
{
    // FIXME: Do we want to get a Nalu authtoken?
    // These can probably stay hard-coded for the foreseeable future
    authtoken_ = "edu.utexas.ices.jmelvin-25970d6f";
    dataset_ = "isotropic1024coarse";
    field_ = "velocity";
 
    // initalize soap service
    soapinit();

    // JHTDB property to close out if error occurs
    turblibSetExitOnError(1);
}

void HITFields::close_JHTDB()
{
   // finalize soap service
    soapdestroy();
}

}  // nalu
}  // sierra
