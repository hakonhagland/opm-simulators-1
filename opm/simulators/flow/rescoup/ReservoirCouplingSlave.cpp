/*
  Copyright 2024 Equinor ASA

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
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <vector>
#include <fmt/format.h>

namespace Opm {
// NOTE: All slave-master communicators have set a custom error handler, which eventually
//   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
//   or MPI_Send() calls below.

template <class Scalar>
ReservoirCouplingSlave<Scalar>::
ReservoirCouplingSlave(
    const Parallel::Communication &comm,
    const Schedule &schedule,
    const SimulatorTimer &timer
) :
    comm_{comm},
    schedule_{schedule},
    timer_{timer}
{
    this->slave_master_comm_ = MPI_COMM_NULL;
    MPI_Comm_get_parent(&this->slave_master_comm_);
    if (this->slave_master_comm_ == MPI_COMM_NULL) {
        OPM_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    }
    // NOTE: By installing a custom error handler for all slave-master communicators, which
    //   eventually will call MPI_Abort(), there is no need to check the return value of any
    //   MPI_Recv() or MPI_Send() calls as errors will be caught by the error handler.
    ReservoirCoupling::setErrhandler(this->slave_master_comm_, /*is_master=*/false);
}

template <class Scalar>
bool
ReservoirCouplingSlave<Scalar>::
isProducer(int report_step_idx) const
{
    return this->schedule_[report_step_idx].getGroup(this->slave_group_order_[report_step_idx]).isProductionGroup();
}

// NOTE: It is not legal for a slave to activate before the master has activated. This problem
//       will be caught by the master when it receives the slave activation date. See:
//       ReservoirCouplingSpawnSlaves::receiveActivationDateFromSlaves_()
template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
maybeActivate(int report_step) {
    if (!this->activated()) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            this->activated_ = true;
            // Send a handshake to the master process to indicate that the slave has activated
            this->sendActivationHandshakeToMasterProcess_();
        }
    }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
maybeReceiveGroupInfoFromMaster() {
    std::size_t size;
    std::vector<char> group_names;
    std::uint8_t start_of_master_report_step_flag;
    if (this->comm_.rank() == 0) {
        auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &start_of_master_report_step_flag,
            /*count=*/1,
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterStartOfReportStep),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info(fmt::format("Received master start of report step flag {} from master process rank 0",
                start_of_master_report_step_flag));
    }
    this->comm_.broadcast(&start_of_master_report_step_flag, /*count=*/1, /*emitter_rank=*/0);
    if (start_of_master_report_step_flag == 1u) {
        this->receiveGroupInfoFromMasterProcess_();
    }
}

template <class Scalar>
double
ReservoirCouplingSlave<Scalar>::
receiveNextTimeStepFromMaster() {
    double timestep;
    if (this->comm_.rank() == 0) {
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &timestep,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextTimeStep),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info(
            fmt::format("Slave rank 0 received next timestep {} from master.", timestep)
        );
    }
    this->comm_.broadcast(&timestep, /*count=*/1, /*emitter_rank=*/0);
    this->logger_.info("Broadcasted slave next time step to all ranks");
    return timestep;
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendAndReceiveInitialData() {
    this->sendActivationDateToMasterProcess_();
    this->sendSimulationStartDateToMasterProcess_();
    this->receiveSlaveNameFromMasterProcess_();
    this->receiveMasterGroupNamesFromMasterProcess_();
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendGroupInfoToMaster(int report_step_idx) const
{
    if (this->comm_.rank() == 0) {
        auto num_slave_groups = this->numSlaveGroups();
        std::vector<std::uint8_t> slave_group_prod_inj_info(num_slave_groups*2);
        for (std::size_t group_idx = 0; group_idx < num_slave_groups; ++group_idx) {
            const auto& group_name = this->slave_group_order_[group_idx];
            const Group& group = this->schedule_.getGroup(group_name, report_step_idx);
            slave_group_prod_inj_info[2*group_idx] = static_cast<std::uint8_t>(group.isInjectionGroup());
            slave_group_prod_inj_info[2*group_idx+1] = static_cast<std::uint8_t>(group.isProductionGroup());
        }
        auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
        MPI_Send(
            slave_group_prod_inj_info.data(),
            /*count=*/slave_group_prod_inj_info.size(),
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveGroupInfo),
            this->slave_master_comm_
        );
    }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendNextReportDateToMasterProcess() const
{
    if (this->comm_.rank() == 0) {
        double elapsed_time = this->timer_.simulationTimeElapsed();
        double current_step_length = this->timer_.currentStepLength();
        // NOTE: This is an offset in seconds from the start date, so it will be 0 if the next report
        //      would be the start date. In general, it should be a positive number.
        double next_report_time_offset = elapsed_time + current_step_length;
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &next_report_time_offset,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextReportDate),
            this->slave_master_comm_
        );
        this->logger_.info("Sent next report date to master process from rank 0");
   }
}


template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendPotentialsToMaster(const std::vector<ReservoirCoupling::Potentials<Scalar>> &potentials) const
{
    // NOTE: The master can determine from the ordering of the potentials in the vector
    //   which slave group for a given slave name the given potentials belong to,
    //   so we do not need to send the slave group names also.
    if (this->comm_.rank() == 0) {
        //auto num_groups = potentials.size();
        auto MPI_POTENTIALS_TYPE = Dune::MPITraits<ReservoirCoupling::Potentials<Scalar>>::getType();
        MPI_Send(
            potentials.data(),
            /*count=*/potentials.size(),
            /*datatype=*/MPI_POTENTIALS_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::Potentials),
            this->slave_master_comm_
        );
        this->logger_.info(
            "Sent potentials to master process from rank 0, slave name: " + this->slave_name_
        );
    }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendStartOfReportStepFlagToMaster(bool start_of_report_step) {
    std::uint8_t start_of_report_step_flag;
    if (this->comm_.rank() == 0) {
        start_of_report_step_flag = start_of_report_step ? 1u : 0u;
        auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &start_of_report_step_flag,
            /*count=*/1,
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveStartOfReportStep),
            this->slave_master_comm_
        );
        this->logger_.info(fmt::format("Sent master start of report step flag {} to master process rank 0",
                start_of_report_step_flag));
    }
    this->comm_.broadcast(&start_of_report_step_flag, /*count=*/1, /*emitter_rank=*/0);
    this->start_of_report_step_flag_ = (start_of_report_step_flag == 1u);
}


// ------------------
// Private methods
// ------------------

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
checkGrupSlavGroupNames_()
{
    // Validate that each slave group name has a corresponding master group name
    bool grup_slav_found = false;
    for (std::size_t report_step = 0; report_step < this->schedule_.size(); ++report_step) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            grup_slav_found = true;
            auto grup_slavs = rescoup.grupSlavs();
            for (const auto& [slave_group_name, grup_slav] : grup_slavs) {
                auto map_iter = this->slave_to_master_group_map_.find(slave_group_name);
                if (map_iter == this->slave_to_master_group_map_.end()) {
                    OPM_THROW(std::runtime_error,
                              "Reservoir coupling: Failed to find master group name for slave group: "
                              + slave_group_name);
                }
                else {
                    const auto& master_group_name = map_iter->second;
                    if (grup_slav.masterGroupName() != master_group_name) {
                        OPM_THROW(std::runtime_error,
                                  "Reservoir coupling: Inconsistent master group name for slave group: "
                                  + slave_group_name);
                    }
                }
            }
        }
    }
    if (!grup_slav_found) {
        OPM_THROW(std::runtime_error, "Reservoir coupling: Failed to find slave group names: "
                  "No GRUPSLAV keyword found in schedule");
    }
}

template <class Scalar>
double
ReservoirCouplingSlave<Scalar>::
getGrupSlavActivationDate_() const
{
    double start_date = this->schedule_.getStartTime();
    for (std::size_t report_step = 0; report_step < this->schedule_.size(); ++report_step) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            return start_date + this->schedule_.seconds(report_step);
        }
    }
    OPM_THROW(std::runtime_error, "Reservoir coupling: Failed to find slave activation time: "
              "No GRUPSLAV keyword found in schedule");
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
receiveGroupInfoFromMasterProcess_() {
    auto num_slave_groups = this->numSlaveGroups();
    std::vector<std::uint8_t> master_group_prod_inj_info_(num_slave_groups*2);
    if (this->comm_.rank() == 0) {
        auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            master_group_prod_inj_info_.data(),
            /*count=*/master_group_prod_inj_info_.size(),
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupInfo),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info(fmt::format(
            "Received master group info for {} groups from master process rank 0", num_slave_groups
        ));
    }
    this->comm_.broadcast(
        master_group_prod_inj_info_.data(), /*count=*/master_group_prod_inj_info_.size(), /*emitter_rank=*/0
    );
    for (std::size_t i = 0; i < num_slave_groups; ++i) {
        this->master_group_is_injector_[this->slave_group_order_[i]] = master_group_prod_inj_info_[2*i];
        this->master_group_is_producer_[this->slave_group_order_[i]] = master_group_prod_inj_info_[2*i+1];
    }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
receiveMasterGroupNamesFromMasterProcess_() {
    std::size_t size;
    std::vector<char> group_names;
    if (this->comm_.rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &size,
            /*count=*/1,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNamesSize),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info("Received master group names size from master process rank 0");
        group_names.resize(size);
        MPI_Recv(
            group_names.data(),
            /*count=*/size,
            /*datatype=*/MPI_CHAR,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNames),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info("Received master group names from master process rank 0");
    }
    this->comm_.broadcast(&size, /*count=*/1, /*emitter_rank=*/0);
    if (this->comm_.rank() != 0) {
        group_names.resize(size);
    }
    this->comm_.broadcast(group_names.data(), /*count=*/size, /*emitter_rank=*/0);
    this->saveMasterGroupNamesAsMapAndEstablishOrder_(group_names);
    this->checkGrupSlavGroupNames_();
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
receiveSlaveNameFromMasterProcess_() {
    std::size_t size;
    std::string slave_name;
    if (this->comm_.rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &size,
            /*count=*/1,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNameSize),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info("Received slave name size from master process rank 0");
        slave_name.resize(size+1); // +1 for the null terminator
        MPI_Recv(
            slave_name.data(),
            /*count=*/size,
            /*datatype=*/MPI_CHAR,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveName),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        slave_name[size] = '\0';  // Add null terminator
        this->logger_.info("Received slave name from master process rank 0");
    }
    this->comm_.broadcast(&size, /*count=*/1, /*emitter_rank=*/0);
    if (this->comm_.rank() != 0) {
        slave_name.resize(size+1); // +1 for the null terminator
    }
    this->comm_.broadcast(slave_name.data(), /*count=*/size+1, /*emitter_rank=*/0);
    this->logger_.info(fmt::format("Received slave name: {}", slave_name));
    this->slave_name_ = slave_name;
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
saveMasterGroupNamesAsMapAndEstablishOrder_(const std::vector<char>& group_names) {
    // Deserialize the group names vector into a map of slavegroup names -> mastergroup names
    // and establish the order the master process sends us the group data. This will enable the
    // master to send group data without the need to send the group names themselves.
    //
    // Call chain: sendAndReceiveInitialData() -> receiveMasterGroupNamesFromMasterProcess_()
    auto total_size = group_names.size();
    std::size_t offset = 0;
    std::size_t idx = 0;
    while (offset < total_size) {
        std::string master_group{group_names.data() + offset};
        offset += master_group.size() + 1;
        assert(offset < total_size);
        std::string slave_group{group_names.data() + offset};
        offset += slave_group.size() + 1;
        this->slave_to_master_group_map_[slave_group] = master_group;
        this->slave_group_order_[idx] = slave_group;
        idx++;
    }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendActivationDateToMasterProcess_() const
{
    if (this->comm_.rank() == 0) {
        // NOTE: The master process will use this date to check that no slave starts before the master
        double activation_date = this->getGrupSlavActivationDate_();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &activation_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveActivationDate),
            this->slave_master_comm_
        );
        this->logger_.info("Sent simulation activation date to master process from rank 0");
   }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendActivationHandshakeToMasterProcess_() const
{
    if (this->comm_.rank() == 0) {
        std::uint8_t activation_handshake = 1u;
        auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &activation_handshake,
            /*count=*/1,
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveActivationHandshake),
            this->slave_master_comm_
        );
        this->logger_.info("Sent simulation activation handshake to master process from rank 0");
    }
    this->comm_.barrier();
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendSimulationStartDateToMasterProcess_() const
{
    if (this->comm_.rank() == 0) {
        // NOTE: The master process needs the s
        double start_date = this->schedule_.getStartTime();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &start_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveSimulationStartDate),
            this->slave_master_comm_
        );
        this->logger_.info("Sent simulation start date to master process from rank 0");
   }
}

// Explicit template instantiations
template class ReservoirCouplingSlave<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingSlave<float>;
#endif

} // namespace Opm
