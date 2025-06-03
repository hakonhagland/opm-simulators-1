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
#include <opm/simulators/wells/GroupTargetCalculator.hpp>
#include <opm/simulators/wells/RescoupTargetCalculator.hpp>

#include <array>
#include <string>
#include <vector>
#include <tuple>

#include <fmt/format.h>

namespace Opm {

// -------------------------------------------------------
// Constructor for the RescoupTargetCalculator class
// -------------------------------------------------------
template <class Scalar>
RescoupTargetCalculator<Scalar>::
RescoupTargetCalculator(
    GuideRateHandler<Scalar>& guide_rate_handler,
    const WellState<Scalar>& well_state,
    const GroupState<Scalar>& group_state,
    const int report_step_idx
)
    : guide_rate_handler_{guide_rate_handler}
    , well_state_{well_state}
    , group_state_{group_state}
    , report_step_idx_{report_step_idx}
    , schedule_{guide_rate_handler.schedule()}
    , summary_state_{guide_rate_handler.summaryState()}
    , deferred_logger_{guide_rate_handler.deferredLogger()}
    , reservoir_coupling_master_{guide_rate_handler.reservoirCouplingMaster()}
    , well_model_{guide_rate_handler.wellModel()}
    , phase_usage_{guide_rate_handler.phaseUsage()}
{
}

template <class Scalar>
void
RescoupTargetCalculator<Scalar>::
calculateMasterGroupTargetsAndSendToSlaves()
{
    // NOTE: Since this object can only be constructed for a master process, we can be
    //   sure that if we are here, we are running as master
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto& comm = rescoup_master.getComm();
    if (comm.rank() == 0) {
        GroupTargetCalculator calculator{
            this->well_model_,
            this->well_state_,
            this->group_state_,
            this->schedule_,
            this->summary_state_,
            this->guide_rate_handler_.phaseUsage(),
            this->guide_rate_handler_.guideRate(),
            this->report_step_idx_,
            this->deferred_logger_
        };
        const auto& slaves = rescoup_master.getSlaveNames();
        static const std::array<Phase, 3> phases = { Phase::WATER, Phase::OIL, Phase::GAS };
        for (const auto& slave_name : slaves) {
            std::vector<InjectionGroupTarget> injection_targets;
            std::vector<ProductionGroupTarget> production_targets;
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_name);
            for (std::size_t idx = 0; idx < master_groups.size(); ++idx) {
                const auto& group_name = master_groups[idx];
                const Group& group = this->schedule_[this->report_step_idx_].groups.get(group_name);
                if (group.isInjectionGroup()) {
                    for (Phase phase : phases) {
                        auto target_info = calculator.groupInjectionTarget(group, phase);
                        if (target_info.has_value()) {
                            injection_targets.push_back(
                                InjectionGroupTarget{
                                    /*group_name_idx=*/idx, target_info->target, target_info->cmode, phase
                                }
                            );
                        }
                    }
                }
                if (group.isProductionGroup()) {
                    auto target_info = calculator.groupProductionTarget(group);
                    if (target_info.has_value()) {
                        production_targets.push_back(
                            ProductionGroupTarget{
                                /*group_name_idx=*/idx, target_info->target, target_info->cmode
                            }
                        );
                    }
                }
            }
            this->sendSlaveGroupTargetsToSlave_(
                rescoup_master, slave_name, injection_targets, production_targets
            );
        }
    }
}

template <class Scalar>
void
RescoupTargetCalculator<Scalar>::
sendSlaveGroupTargetsToSlave_(
    const ReservoirCouplingMaster<Scalar>& rescoup_master,
    const std::string& slave_name,
    const std::vector<InjectionGroupTarget>& injection_targets,
    const std::vector<ProductionGroupTarget>& production_targets
) const
{
    auto num_injection_targets = injection_targets.size();
    auto num_production_targets = production_targets.size();
    // First, send the number of targets such that the slave can know if it can expect none
    // or more targets.
    rescoup_master.sendNumGroupTargetsToSlave(slave_name, num_injection_targets, num_production_targets);
    if (num_injection_targets > 0) {
        rescoup_master.sendInjectionTargetsToSlave(slave_name, injection_targets);
    }
    if (num_production_targets > 0) {
        rescoup_master.sendProductionTargetsToSlave(slave_name, production_targets);
    }
}

template class RescoupTargetCalculator<double>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupTargetCalculator<float>;
#endif

}// namespace Opm
