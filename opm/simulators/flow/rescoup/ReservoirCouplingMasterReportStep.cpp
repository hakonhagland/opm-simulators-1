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

#include <config.h>
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <vector>
#include <fmt/format.h>

namespace Opm {

template <class Scalar>
ReservoirCouplingMasterReportStep<Scalar>::
ReservoirCouplingMasterReportStep(
    ReservoirCouplingMaster<Scalar> &master
) :
    master_{master}
{
}

// ------------------
// Public methods
// ------------------

template <class Scalar>
const ReservoirCoupling::Potentials<Scalar>&
ReservoirCouplingMasterReportStep<Scalar>::
getSlaveGroupPotentials(const std::string &master_group_name)
{
    auto it = this->master_group_slave_names_.find(master_group_name);
    if (it != this->master_group_slave_names_.end()) {
        auto& slave_name = it->second;
        auto pot_idx = this->getMasterGroupCanonicalIdx(slave_name, master_group_name);
        return this->slave_group_potentials_[slave_name][pot_idx];
    }
    else {
        OPM_THROW(
            std::runtime_error,
            fmt::format(
                "Master group name {} not found in master-to-slave-group-name mapping",
                master_group_name
            )
        );
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
maybeReceiveGroupInfoFromSlaves()
{
    auto num_slaves = this->numSlaves();
    this->logger().info("Receiving group info from slave processes");
    for (unsigned int slave_idx = 0; slave_idx < num_slaves; slave_idx++) {
        if (!this->slaveIsActivated(slave_idx)) {
            continue;
        }
        if (this->comm().rank() == 0) {
            std::uint8_t start_of_report_step_flag;
            auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
            // NOTE: See comment about error handling at the top of ReservoirCouplingMaster.cpp
            MPI_Recv(
                &start_of_report_step_flag,
                /*count=*/1,
                /*datatype=*/MPI_UINT8_T_TYPE,
                /*source_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveStartOfReportStep),
                this->getSlaveComm(slave_idx),
                MPI_STATUS_IGNORE
            );
            this->logger().info(fmt::format(
                "Received start of report step flag {} from slave process with name: {}",
                start_of_report_step_flag, this->slaveName(slave_idx)
            ));
        }
        this->comm().broadcast(&start_of_report_step_flag, /*count=*/1, /*emitter_rank=*/0);
        if (start_of_report_step_flag == 1u) {
            this->receiveGroupInfoFromSlave_(slave_idx);
        }
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
receivePotentialsFromSlaves()
{
    auto num_slaves = this->numSlavesStarted();
    this->logger_.info("Receiving potentials from slave processes");
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto num_slave_groups = this->numSlaveGroups(i);
        assert( num_slave_groups > 0 );
        std::vector<ReservoirCoupling::Potentials<Scalar>> potentials(num_slave_groups);
        if (this->comm_.rank() == 0) {
            if (this->slaveIsActivated(i)) {
                // NOTE: See comment about error handling at the top of ReservoirCouplingMaster.cpp
                auto MPI_POTENTIALS_TYPE = Dune::MPITraits<ReservoirCoupling::Potentials<Scalar>>::getType();
                MPI_Recv(
                    potentials.data(),
                    /*count=*/num_slave_groups,
                    /*datatype=*/MPI_POTENTIALS_TYPE,
                    /*source_rank=*/0,
                    /*tag=*/static_cast<int>(MessageTag::Potentials),
                    this->getSlaveComm(i),
                    MPI_STATUS_IGNORE
                );
                this->logger_.info(
                    fmt::format(
                        "Received potentials from slave process with name: {}. "
                        "Number of slave groups: {}", this->slave_names_[i], num_slave_groups
                    )
                );
            }
            else {
                this->logger_.info(fmt::format(
                    "Slave {} has not activated yet, skipping receiving potentials from slave",
                        this->slave_names_[i]
                    )
                );
                potentials.assign(num_slave_groups, ReservoirCoupling::Potentials<Scalar>{}); // Set to zero potentials
            }
        }
        // NOTE: The dune broadcast() below will do something like:
        //    MPI_Bcast(inout,len,MPITraits<Potentials>::getType(),root,communicator)
        //  so it should use the custom Potentials MPI type that we defined in
        //  ReservoirCouplingMpiTraits.hpp
        this->comm_.broadcast(potentials.data(), /*count=*/num_slave_groups, /*emitter_rank=*/0);
        this->slave_group_potentials_[this->slave_names_[i]] = potentials;
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
sendGroupInfoToSlaves(int report_step_idx)
{
    if (this->comm().rank() == 0) {
        for (unsigned int slave_idx = 0; slave_idx < this->numSlaves(); slave_idx++) {
            if (!this->slaveIsActivated(slave_idx)) {
                this->logger().info(fmt::format(
                    "Slave {} has not activated yet, skipping sending group info to slave",
                    this->slaveName(slave_idx)
                ));
                continue;
            }
            this->sendMasterGroupInjectorProducerInfoToSlave_(slave_idx, report_step_idx);
        }
   }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
setReportStepIdx(int report_step_idx)
{
    this->report_step_idx_ = report_step_idx;
}


// ------------------
// Private methods
// ------------------

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
receiveGroupInfoFromSlave_(unsigned int slave_idx)
{
    auto num_slave_groups = this->numSlaveGroups(slave_idx);
    const auto& slave_name = this->slaveName(slave_idx);
    assert( num_slave_groups > 0 );
    std::vector<std::uint8_t> group_info(2*num_slave_groups);
    // NOTE: See comment about error handling at the top of ReservoirCouplingMaster.cpp
    if (this->comm().rank() == 0) {
        MPI_Recv(
            group_info.data(),
            /*count=*/2*num_slave_groups,
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveGroupInfo),
            this->getSlaveComm(slave_idx),
            MPI_STATUS_IGNORE
        );
        this->logger().info(fmt::format(
            "Received group info from slave process with name: {}", slave_name
        ));
    }
    this->comm().broadcast(group_info.data(), /*count=*/2*num_slave_groups, /*emitter_rank=*/0);
    // Resize the vectors to the number of slave groups first time (later this should be a no-op)
    this->slave_group_is_producer_[slave_name].resize(num_slave_groups);
    this->slave_group_is_injector_[slave_name].resize(num_slave_groups);
    // Set the group info
    for (std::size_t group_idx = 0; group_idx < num_slave_groups; ++group_idx) {
        this->slave_group_is_injector_[slave_name][group_idx] = group_info[2*group_idx];
        this->slave_group_is_producer_[slave_name][group_idx] = group_info[2*group_idx+1];
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
sendMasterGroupInjectorProducerInfoToSlave_(unsigned int slave_idx, int report_step_idx)
{
    // NOTE: This vector have an order which is known by the slave process, so by sending the group info
    //       in the same order, the slave process can reconstruct the master group (slave group) names
    //       avoiding any unnecessary communication of group names.
    const auto& master_groups = this->getMasterGroupNamesForSlave(slave_idx);
    std::size_t num_master_groups = master_groups.size();
    std::vector<std::uint8_t> group_info(2*num_master_groups);
    for (std::size_t group_idx = 0; group_idx < num_master_groups; ++group_idx) {
        const auto& group_name = master_groups[group_idx];
        const Group& group = this->schedule().getGroup(group_name, report_step_idx);
        group_info[2*group_idx] = static_cast<std::uint8_t>(group.isInjectionGroup());
        group_info[2*group_idx + 1] = static_cast<std::uint8_t>(group.isProductionGroup());
    }
    auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
    // NOTE: See comment about error handling at the top of ReservoirCouplingMaster.cpp
    MPI_Send(
        group_info.data(),
        /*count=*/group_info.size(),
        /*datatype=*/MPI_UINT8_T_TYPE,
        /*dest_rank=*/0,
        /*tag=*/static_cast<int>(MessageTag::MasterGroupInfo),
        this->getSlaveComm(slave_idx)
    );
    this->logger().info(fmt::format(
        "Sent master group info from master process rank 0 to slave process "
        "rank 0 with name: {}", this->slaveName(slave_idx))
    );
}

// Explicit instantiations
template class ReservoirCouplingMasterReportStep<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingMasterReportStep<float>;
#endif

} // namespace Opm
