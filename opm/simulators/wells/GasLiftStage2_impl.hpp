/*
  Copyright 2020 Equinor ASA.

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

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>

#include <optional>
#include <string>

#include <fmt/format.h>

template<typename TypeTag>
Opm::GasLiftStage2<TypeTag>::
GasLiftStage2(
    const BlackoilWellModel &well_model,
    const Simulator &ebos_simulator,
    DeferredLogger &deferred_logger,
    const WellState &well_state,
    GasLiftWells &glift_wells
) :
    deferred_logger_{deferred_logger},
    ebos_simulator_{ebos_simulator},
    well_model_{well_model},
    well_state_{well_state},
    stage1_wells_{glift_wells},
    report_step_idx_{ebos_simulator_.episodeIndex()},
    summary_state_{ebos_simulator_.vanguard().summaryState()},
    comm_{ebos_simulator_.vanguard().grid().comm()},
    schedule_{ebos_simulator.vanguard().schedule()},
    phase_usage_{well_model_.phaseUsage()}
{ }

/****************************************
 * Methods in alphabetical order
 ****************************************/

template<typename TypeTag>
void
Opm::GasLiftStage2<TypeTag>::
displayDebugMessage_(const std::string &msg, const std::string &group_name)
{
    const std::string message = fmt::format(
        "  GLIFT2 (DEBUG) : Group {} : {}", group_name, msg);
    this->deferred_logger_.info(message);
}

// Find all subordinate wells of a given group.
//
// NOTE: A group can either contain wells or groups, not both.
//   If it contains groups, we have to traverse those recursively to find the wells.
//
template<typename TypeTag>
std::vector<Opm::GasLiftSingleWell<TypeTag> *>
Opm::GasLiftStage2<TypeTag>::
getGroupGliftWells_(const Opm::Group &group)
{
    std::vector<GasLiftSingleWell *> wells;
    getGroupGliftWellsRecursive_(group, wells);
    return wells;

}

template<typename TypeTag>
void
Opm::GasLiftStage2<TypeTag>::
getGroupGliftWellsRecursive_(const Opm::Group &group,
         std::vector<GasLiftSingleWell*> &wells)
{
    for (const std::string& group_name : group.groups()) {
        const Group& sub_group = this->schedule_.getGroup(
            group_name, this->report_step_idx_);
        getGroupGliftWellsRecursive_(sub_group, wells);
    }
    for (const std::string& well_name : group.wells()) {
        if (this->stage1_wells_.count(well_name) == 1) {
            GasLiftSingleWell *well_ptr = this->stage1_wells_.at(well_name).get();
            wells.push_back(well_ptr);
        }
    }
}

template<typename TypeTag>
void
Opm::GasLiftStage2<TypeTag>::
optimizeGroup_(const Opm::Group &group)
{
    const GasLiftOpt& glo = this->schedule_.glo(this->report_step_idx_);
    for (const std::string& group_name : group.groups()) {
        const Group& sub_group = this->schedule_.getGroup(
            group_name, this->report_step_idx_);
        optimizeGroup_(sub_group);
    }
    try {
        /*const auto &gl_group =*/ glo.group(group.name());
    }
    catch (std::out_of_range &e) {
        displayDebugMessage_("no gaslift info available", group.name());
        return;
    }
    const auto &gl_group = glo.group(group.name());
    const auto &max_glift = gl_group.max_lift_gas();
    if (group.has_control(Group::ProductionCMode::ORAT) || max_glift) {
        //const auto controls = group.productionControls(this->summary_state_);
        //const auto &rates = WellGroupHelpers::getProductionGroupRateVector(
        //    this->well_state_, this->phase_usage_, group.name());
        //if ((controls.oil_target < rates.oil_rat) || max_glift ) {.. }
        displayDebugMessage_("optimizing", group.name());
        for (auto well : getGroupGliftWells_(group)) {
            const std::string msg = fmt::format("well: {}", well->name());
            displayDebugMessage_(msg, group.name());

                /*
                const auto& it = this->well_state_.wellMap().find(well_name);
                if (it == end) continue;
                int well_index = it->second[0];
                const auto& well_ecl = this->schedule_.getWell(
                    well_name, this->report_step_idx_);
                //We don't have to check if the stage1_wells are producers
                //  or if they have THP control..
                if (this->stage1_wells.count(well_name)) {
                if (well_ecl.isProducer()) {
                }

            if (wellEcl.getStatus() == Well::Status::SHUT)
                continue;

            double factor = wellEcl.getEfficiencyFactor();

            }
                */
        }
    }
    else {
        displayDebugMessage_("skipping", group.name());
    }
}

template<typename TypeTag>
void
Opm::GasLiftStage2<TypeTag>::
runOptimize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);

    optimizeGroup_(group);

}
