/*
  Copyright 2024 Equinor AS

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

#include <config.h>
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <vector>

namespace Opm {

ReservoirCouplingSlave::ReservoirCouplingSlave(
    const Parallel::Communication &comm,
    const Schedule &schedule,
    const SimulatorTimer &timer
) :
    comm_{comm},
    schedule_{schedule},
    timer_{timer}
{
    this->slave_master_comm_ = MPI_Comm_Ptr(new MPI_Comm(MPI_COMM_NULL));
    MPI_Comm_get_parent(this->slave_master_comm_.get());
    if (*(this->slave_master_comm_) == MPI_COMM_NULL) {
        OPM_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    }
}

void ReservoirCouplingSlave::sendNextReportDateToMasterProcess() {
    this->comm_.barrier();  // For debugging purposes
    OpmLog::info("xxx6: comm.size = " + std::to_string(this->comm_.size()));
    OpmLog::info("xxx5: Rank " + std::to_string(this->comm_.rank()) + " reached the start of the method.");
    if (this->comm_.rank() == 0) {
        double elapsed_time = this->timer_.simulationTimeElapsed();
        double current_step_length = this->timer_.currentStepLength();
        double next_report_date = elapsed_time + current_step_length;
        OpmLog::info("xxx2: Sending next report date..");
        MPI_Send(
            &next_report_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextReportDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent next report date to master process from rank 0");
   }
   OpmLog::info("xxx4: Rank " + std::to_string(this->comm_.rank()) + " reached the barrier.");
   this->comm_.barrier();  // For debugging purposes, to ensure that the master process receives the message
   OpmLog::info("xxx3: Sent next report date to master process from rank 0");
}

void ReservoirCouplingSlave::sendSimulationStartDateToMasterProcess() {
    OpmLog::info("xxx7: Rank " + std::to_string(this->comm_.rank()) + " reached the start of the method.");
    std::cout << "xxx8: Rank " << this->comm_.rank() << " reached the start of the method." << std::endl;
    if (this->comm_.rank() == 0) {
        // Ensure that std::time_t is of type long since we are sending it over MPI with MPI_LONG
        static_assert(std::is_same<std::time_t, long>::value, "std::time_t is not of type long");
        std::time_t start_date = this->schedule_.getStartTime();
        MPI_Send(
            &start_date,
            /*count=*/1,
            /*datatype=*/MPI_LONG,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveSimulationStartDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent simulation start date to master process from rank 0");
   }
}

} // namespace Opm
