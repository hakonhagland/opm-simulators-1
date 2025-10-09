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
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/simulators/wells/rescoup/RescoupReceiveSlaveGroupData.hpp>

#include <array>
#include <string>
#include <vector>
#include <tuple>

#include <fmt/format.h>

namespace Opm {

// -------------------------------------------------------
// Constructor for the RescoupTargetCalculator class
// -------------------------------------------------------
template <class Scalar, class IndexTraits>
RescoupReceiveSlaveGroupData<Scalar, IndexTraits>::
RescoupReceiveSlaveGroupData(
    ReservoirCouplingMaster<Scalar>& reservoir_coupling_master,
    const Schedule& schedule,
    const GroupState<Scalar>& group_state,
    const PhaseUsageInfo<IndexTraits>& phase_usage,
    const int report_step_idx
)
    : reservoir_coupling_master_{reservoir_coupling_master}
    , schedule_{schedule}
    , group_state_{group_state}
    , phase_usage_{phase_usage}
    , report_step_idx_{report_step_idx}
{
}

template <class Scalar, class IndexTraits>
void
RescoupReceiveSlaveGroupData<Scalar, IndexTraits>::
receiveSlaveGroupData()
{
    auto& rescoup_master = this->reservoir_coupling_master_;
    // Compute size of production data based on the number of slave groups that are producers
    auto num_slave_groups = rescoup_master.numSlaveGroups();
    std::size_t num_production_data = 0;
    for (std::size_t slave_idx = 0; slave_idx < num_slave_groups; ++slave_idx) {
        std::size_t master_group_idx = 0;
        for (const auto& master_group_name : rescoup_master.getMasterGroupNamesForSlave(slave_idx)) {
            const Group& group = this->schedule_.getGroup(master_group_name, this->report_step_idx_);
            if (group.isProducer(this->report_step_idx_)) {
                if (rescoup_master.slaveGroupIsProducer(slave_idx, master_group_idx)) {
                    num_production_data++;
                }
                master_group_idx++;
            }
        }
    this->receiveSlaveGroupProductionData_();
        this->receiveSlaveGroupInjectionData_();
    }
}


// -------------------------------------------------------
// Private methods for class RescoupReceiveSlaveGroupData
// -------------------------------------------------------
template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::Potentials
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupPotentials_(std::size_t group_idx)
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    Potentials potentials;
    assert(group.isProducer(this->report_step_idx_));
    const auto& gr_pot = this->group_state_.get_production_group_potential(group_name);
    potentials[Potentials::Phase::Oil] = gr_pot.oil_rate;
    potentials[Potentials::Phase::Gas] = gr_pot.gas_rate;
    potentials[Potentials::Phase::Water] = gr_pot.water_rate;
    return potentials;
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::ProductionRates
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupProductionRates_(std::size_t group_idx)
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    assert(group.isProducer(this->report_step_idx_));
    GuideRate::RateVector production_rates = WellGroupHelpersType::getProductionGroupRateVector(
        this->group_state_,
        this->phase_usage_,
        group_name
    );
    // NOTE: GuideRate::RateVector is a vector of doubles, so we need to convert it to Scalar
    // TODO: Fix GuideRate::RateVector to be a vector of Scalars instead of doubles
    return ProductionRates{production_rates};
}


template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::SlaveGroupProductionData
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupProductionData_(std::size_t group_idx)
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    SlaveGroupProductionData production_data;
    // Potentials are used to calculate guide rates for master production groups
    production_data.potentials = this->collectSlaveGroupPotentials_(group_idx);
    // Production rates are used to transform guiderate targets for master production groups
    //   from one phase to another
    production_data.production_rates = this->collectSlaveGroupProductionRates_(group_idx);
    return production_data;
}

template <class Scalar, class IndexTraits>
void
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
sendSlaveGroupProductionDataToMaster_()
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    auto num_slave_groups = rescoup_slave.numSlaveGroups();
    std::vector<SlaveGroupProductionData> production_data;
    for (std::size_t group_idx = 0; group_idx < num_slave_groups; ++group_idx) {
        if (rescoup_slave.masterGroupIsProducer(group_idx)) {
            const auto& group_name = rescoup_slave.slaveIdxToGroupName(group_idx);
            const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
            if (group.isProducer(this->report_step_idx_)) {
                production_data.emplace_back(this->collectSlaveGroupProductionData_(group_idx));
            }
        }
    }
    rescoup_slave.sendSlaveGroupProductionDataToMaster(production_data);
}

template <class Scalar, class IndexTraits>
void
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
sendSlaveGroupInjectionDataToMaster_()
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    auto num_slave_groups = rescoup_slave.numSlaveGroups();
    std::vector<SlaveGroupInjectionData> injection_data;
    for (std::size_t group_idx = 0; group_idx < num_slave_groups; ++group_idx) {
        if (rescoup_slave.masterGroupIsInjector(group_idx)) {
            const auto& group_name = rescoup_slave.slaveIdxToGroupName(group_idx);
            const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
            if (group.isInjector(this->report_step_idx_)) {
                injection_data.emplace_back(this->collectSlaveGroupInjectionData_(group_idx));
            }
        }
    }
    rescoup_slave.sendSlaveGroupInjectionDataToMaster(injection_data);
}

} // namespace Opm

template class RescoupSendSlaveGroupData<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupSendSlaveGroupData<float, BlackOilDefaultFluidSystemIndices>;
#endif
