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
#include <opm/simulators/wells/rescoup/RescoupSendSlaveGroupData.hpp>

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
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
RescoupSendSlaveGroupData(
    ReservoirCouplingSlave<Scalar>& reservoir_coupling_slave,
    const Schedule& schedule,
    const GroupState<Scalar>& group_state,
    const PhaseUsageInfo<IndexTraits>& phase_usage,
    const int report_step_idx
)
    : reservoir_coupling_slave_{reservoir_coupling_slave}
    , schedule_{schedule}
    , group_state_{group_state}
    , phase_usage_{phase_usage}
    , report_step_idx_{report_step_idx}
{
}

// ------------------
// Public methods
// ------------------

// None, some, or all of the following data should be sent to the master process at the beginning of the
// timestep depending on the situation:
// 1) *Production group potentials*: If both the master group and slave groups are production groups (and
//   possibly also injection groups), and the master group has GCONPROD item 10 set to "FORM", "POTN", or "0".
//   In this case, the master group needs the slave group potentials to calculate the guide rates.
//   TODO: Investigate if this could also apply to a parent of the master (and not just the master group itself)
//    For example, if a parent group of the master have GCONPROD item 10 set to "FORM", but the master group
//    itself has no guidrate.
//
// 2) *Production group surface rates*: If both the master group and slave groups are production
//   groups (and possibly also injection groups), and the master group has GCONPROD item 10 set to
//   OIL, GAS, WATER, or LIQ and the phase under group control is different from the phase in GCONPROD item 10.
//   In this case, the master group needs the slave group production rates to transform guide rate targets to
//   the phase under group control. NOTE: This should also apply to parents of the master group.
//
// 3) *Production group surface rate for phase P for reinjection*: If the master group is both an injector
//   and a producer group, and the slave group is a producer group (and possibly also an injection group), and
//   the master group has GCONINJE item 3 set to REIN for a phase P (set in GCONINJE item 2). In this case,
//   the master group needs the slave group reinjection surface rates to calculate the injection target for
//   phase P. NOTE1: If the phase P is GAS, the reinjection surface rate sent from the slave group also
//   includes the gas import and consumption rates. NOTE2: Also applies to a parent of the master group.
//
// 4) *Production group reservoir rates*: a) If the master group is both an injector and a producer group, and
//   the slave group is a producer group (and possibly also an injection group) and the master group has
//   GCONINJE item 10 set to NETV or VOID. In this case, the master group needs the slave group reservoir rates to
//   calculate the guide rates. NOTE: Also applies to parents of the master group.
//
// 5) If the master group or any of its parents is a pressure maintenance group. It may need
//   a) slave group reservoir production rates (GPMAINT item 2 = PROD)
//   b) slave group reservoir injection rates (GPMAINT item 2 = OINJ, WINJ, or GINJ)
//   c) slave group surface injection reates (GPMAINT item 2 = OINS, WINS, or GINS)
//   
//   is both an injector and a producer group, and the slave group is a producer group (and possibly also
//   an injection group) and the master group has GPMAINT
//
//   the master group has GCONINJ item 10 set to RATE, VOID, NETV, or RESV. In this case, the master group
//   needs the slave group injection surface rates to calculate the guide rates.
// - Injection group injection reservoir rates: If both the master group and slave groups are injection groups and
//   the master group has GCONINJ item 10 set to RESV. In this case, the master group needs the slave group
//   injection reservoir rates to calculate the guide rates.
template <class Scalar, class IndexTraits>
void
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
sendSlaveGroupDataToMaster()
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    if (rescoup_slave.comm().rank() == 0) {
        this->sendSlaveGroupProductionDataToMaster_();
        this->sendSlaveGroupInjectionDataToMaster_();
    }
}


// -------------------------------------------------------
// Private methods for class RescoupSendSlaveGroupData
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
            const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
            const Group& slave_group = this->schedule_.getGroup(group_name, this->report_step_idx_);
            if (slave_group.isProducer(this->report_step_idx_)) {
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
