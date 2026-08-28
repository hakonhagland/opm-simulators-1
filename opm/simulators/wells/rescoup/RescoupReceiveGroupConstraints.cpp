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
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/simulators/wells/rescoup/RescoupReceiveGroupConstraints.hpp>

#include <fmt/format.h>

namespace Opm {

template <class Scalar, class IndexTraits>
RescoupReceiveGroupConstraints<Scalar, IndexTraits>::
RescoupReceiveGroupConstraints(
    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
    GroupStateHelper<Scalar, IndexTraits>& group_state_helper
)
    : guide_rate_handler_{guide_rate_handler}
    , group_state_helper_{group_state_helper}
    , reservoir_coupling_slave_{guide_rate_handler.reservoirCouplingSlave()}
{
}

template <class Scalar, class IndexTraits>
void
RescoupReceiveGroupConstraints<Scalar, IndexTraits>::
receiveGroupConstraintsFromMaster()
{
    // NOTE: All ranks must call these functions because they contain broadcasts.
    //   The MPI_Recv parts inside the functions have their own rank 0 checks.
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    auto [num_inj_targets, num_prod_constraints] = rescoup_slave.receiveNumGroupConstraintsFromMaster();
    auto& rc_summary_state = rescoup_slave.reservoirCouplingSummaryState();
    this->group_state_helper_.setReservoirCouplingSummaryState(&rc_summary_state);
    if (num_inj_targets > 0) {
        rescoup_slave.receiveInjectionGroupTargetsFromMaster(num_inj_targets);
    }
    for (std::size_t group_idx = 0; group_idx < rescoup_slave.numSlaveGroups(); ++group_idx) {
        const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
        for (const auto phase : {Phase::WATER, Phase::OIL, Phase::GAS}) {
            if (!rescoup_slave.hasMasterInjectionTarget(group_name, phase)) {
                continue;
            }

            const auto target = rescoup_slave.masterInjectionTarget(group_name, phase).first;
            const auto summary_keyword = phase == Phase::WATER ? "GWIRT" :
                                         phase == Phase::OIL   ? "GOIRT" : "GGIRT";
            const auto measure = phase == Phase::GAS
                ? UnitSystem::measure::gas_surface_rate
                : UnitSystem::measure::liquid_surface_rate;
            const auto target_raw = this->group_state_helper_.schedule().getUnits()
                .from_si(measure, target);
            rc_summary_state.setGroupValue(group_name, summary_keyword, target_raw);
        }
    }
    if (num_prod_constraints > 0) {
        rescoup_slave.receiveProductionGroupConstraintsFromMaster(num_prod_constraints);
    }
}

template class RescoupReceiveGroupConstraints<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupReceiveGroupConstraints<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
