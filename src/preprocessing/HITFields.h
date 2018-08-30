//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#ifndef HITFIELDS_H
#define HITFIELDS_H

#include "PreProcessingTask.h"
#include <stdio.h>
#include <math.h>
#include "JHTDB/turblib.h"

#define PI 3.141592653589793

namespace sierra {
namespace nalu {

class HITFields: public PreProcessingTask
{
public:
    /**
     * \param mesh A sierra::nalu::CFDMesh instance
     * \param node The YAML::Node containing inputs for this task
     */
    HITFields(CFDMesh&, const YAML::Node&);

    virtual ~HITFields() {}

    //! Declare velocity and temperature fields and register them for output
    void initialize();

    //! Initialize the velocity and/or temperature fields by linear interpolation
    void run();

private:
    HITFields() = delete;
    HITFields(const HITFields&) = delete;

    size_t fh_;

    //! Parse the YAML file and initialize parameters
    void load(const YAML::Node&);

    stk::mesh::PartVector fluid_parts_;

    std::vector<double> mean_vel_{0.0, 0.0, 0.0};

    double ti_;
    double urms_;

    std::vector<double> physLen_{0.0, 0.0, 0.0};
    std::vector<int> numPts_{0, 0, 0};

    double minx_;
    double t0_;
    double tf_;
    double filtWidth_;
    int nx_;
    int ny_;
    int nz_;
    double uSweep_;
    double Lx_;
    double Ly_;
    double Lz_;
    int numSteps_;
    double dtPhysical_;

    std::string authtoken_;
    std::string dataset_;
    std::string field_;

    void init_JHTDB();
    void close_JHTDB();

};

}  // nalu
}  // sierra


#endif /* HITFIELDS_H */
