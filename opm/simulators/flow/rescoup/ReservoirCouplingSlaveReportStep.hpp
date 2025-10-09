/*
  Copyright 2025 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_RESERVOIR_COUPLING_SLAVE_REPORT_STEP_HPP
#define OPM_RESERVOIR_COUPLING_SLAVE_REPORT_STEP_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

template <class Scalar>
class ReservoirCouplingSlaveReportStep {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;

    ReservoirCouplingSlaveReportStep(
        const Parallel::Communication &comm, const Schedule &schedule, const SimulatorTimer &timer
    );
};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_SLAVE_REPORT_STEP_HPP
