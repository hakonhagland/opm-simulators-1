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
GasLiftStage2<TypeTag>::
GasLiftStage2(
    const BlackoilWellModel &well_model,
    const Simulator &ebos_simulator,
    DeferredLogger &deferred_logger,
    const WellState &well_state,
    GLiftProdWells &prod_wells,
    GLiftOptWells &glift_wells
) :
    deferred_logger_{deferred_logger},
    ebos_simulator_{ebos_simulator},
    well_model_{well_model},
    well_state_{well_state},
    prod_wells_{prod_wells},
    stage1_wells_{glift_wells},
    report_step_idx_{ebos_simulator_.episodeIndex()},
    summary_state_{ebos_simulator_.vanguard().summaryState()},
    comm_{ebos_simulator_.vanguard().grid().comm()},
    schedule_{ebos_simulator.vanguard().schedule()},
    phase_usage_{well_model_.phaseUsage()},
    glo_{schedule_.glo(report_step_idx_)},
    debug_{true}
{ }

/****************************************
 * Methods in alphabetical order
 ****************************************/

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
deleteGrad_(GradMap &map, const std::string &name)
{
    auto value = map.at(name);
    map.erase(name);
    return value;
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
deleteDecGradItem_(const std::string &name)
{
    return deleteGrad_(this->dec_grads_, name);
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
deleteIncGradItem_(const std::string &name)
{
    return deleteGrad_(this->inc_grads_, name);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayWarning_(const std::string &msg, const std::string &group_name)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION (STAGE2), GROUP: {} : {}", group_name, msg);
    this->deferred_logger_.warning("WARNING", message);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayDebugMessage_(const std::string &msg)
{
    if (this->debug_) {
        const std::string message = fmt::format(
            "  GLIFT2 (DEBUG) : {}", msg);
        this->deferred_logger_.info(message);
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayDebugMessage_(const std::string &msg, const std::string &group_name)
{
    if (this->debug_) {
        const std::string message = fmt::format(
            "Group {} : {}", group_name, msg);
        displayDebugMessage_(message);
    }
}

// Find all subordinate wells of a given group.
//
// NOTE: A group can either contain wells or groups, not both.
//   If it contains groups, we have to traverse those recursively to find the wells.
//
template<typename TypeTag>
std::vector<GasLiftSingleWell<TypeTag> *>
GasLiftStage2<TypeTag>::
getGroupGliftWells_(const Opm::Group &group)
{
    std::vector<GasLiftSingleWell *> wells;
    getGroupGliftWellsRecursive_(group, wells);
    return wells;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
getGroupGliftWellsRecursive_(const Opm::Group &group,
         std::vector<GasLiftSingleWell *> &wells)
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
GasLiftStage2<TypeTag>::
optimizeGroup_(const Opm::Group &group)
{
    const auto &gl_group = this->glo_.group(group.name());
    const auto &max_glift = gl_group.max_lift_gas();
    if (group.has_control(Group::ProductionCMode::ORAT) || max_glift) {
        displayDebugMessage_("optimizing", group.name());
        auto wells = getGroupGliftWells_(group);
        OptimizeState state {*this, group};
        if (!state.checkAtLeastTwoWells(wells)) {
            return;
        }
        assert(wells.size() >= 2);
        std::vector<GradPair> inc_grads;
        std::vector<GradPair> dec_grads;
        inc_grads.reserve(wells.size());
        dec_grads.reserve(wells.size());
        state.calculateEcoGradients(wells, inc_grads, dec_grads);
        bool stop_iteration = false;
        while (!stop_iteration && (state.it++ <= this->max_iterations_)) {
            state.debugShowIterationInfo();
            auto [min_dec_grad, max_inc_grad]
                = state.getEcoGradients(inc_grads, dec_grads);
            if (min_dec_grad) {
                assert( max_inc_grad );
                if ((*max_inc_grad)->second > (*min_dec_grad)->second) {
                    state.redistributeALQ(*min_dec_grad, *max_inc_grad);
                    state.recalculateGradients(
                        inc_grads, dec_grads, *min_dec_grad, *max_inc_grad);
                    continue;
                }
            }
            stop_iteration = true;
        }
        if (state.it > this->max_iterations_) {
            const std::string msg = fmt::format("Max iterations {} exceeded.",
                this->max_iterations_);
            displayWarning_(msg, group.name());
        }
    }
    else {
        displayDebugMessage_("skipping", group.name());
    }

}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
optimizeGroupsRecursive_(const Opm::Group &group)
{
    for (const std::string& group_name : group.groups()) {
        const Group& sub_group = this->schedule_.getGroup(
            group_name, this->report_step_idx_);
        optimizeGroupsRecursive_(sub_group);
    }
    try {
        /*const auto &gl_group =*/ this->glo_.group(group.name());
    }
    catch (std::out_of_range &e) {
        //displayDebugMessage_("no gaslift info available", group.name());
        return;
    }
    optimizeGroup_(group);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
runOptimize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);

    optimizeGroupsRecursive_(group);

}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
saveGrad_(GradMap &map, const std::string &name, GradInfo &grad)
{
    if (auto it = map.find(name); it == map.end()) {
        auto [map_it, success] = map.emplace(name, grad);
        assert(success);
    }
    else {
        it->second = grad;
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
saveDecGrad_(const std::string &name, GradInfo &grad)
{
    saveGrad_(this->dec_grads_, name, grad);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
saveIncGrad_(const std::string &name, GradInfo &grad)
{
    saveGrad_(this->inc_grads_, name, grad);
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
updateGrad_(GradMap &map, const std::string &name, GradInfo &grad)
{
    auto old_value = map.at(name);
    map[name] = grad;
    return old_value;
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
updateDecGrad_(const std::string &name, GradInfo &grad)
{
    return updateGrad_(this->dec_grads_, name, grad);
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
updateIncGrad_(const std::string &name, GradInfo &grad)
{
    return updateGrad_(this->inc_grads_, name, grad);
}

/****************************************
 * Methods declared in OptimizeState
 ****************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
calculateEcoGradients(std::vector<GasLiftSingleWell *> &wells,
           std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    assert(wells.size() >= 2);
    for (auto well_ptr : wells) {
        const auto &gs_well = *well_ptr;  // gs = GasLiftSingleWell
        const auto &name = gs_well.name();
        auto inc_grad = calcIncOrDecGrad_(name, gs_well, /*increase=*/true);
        if (inc_grad) {
            inc_grads.emplace_back(std::make_pair(name, inc_grad->grad));
            this->parent.saveIncGrad_(name, *inc_grad);
        }
        auto dec_grad = calcIncOrDecGrad_(name, gs_well, /*increase=*/false);
        if (dec_grad) {
            dec_grads.emplace_back(std::make_pair(name, dec_grad->grad));
            this->parent.saveDecGrad_(name, *dec_grad);
        }
    }
}


template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::OptimizeState::
checkAtLeastTwoWells(std::vector<GasLiftSingleWell *> &wells)
{
    if (wells.size() < 2) {
        const std::string msg = fmt::format(
            "skipping: too few wells ({}) to optimize.", wells.size());
        displayDebugMessage_(msg);
        return false;
    }
    return true;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
debugShowIterationInfo()
{
    const std::string msg = fmt::format("iteration {}", this->it);
    displayDebugMessage_(msg);
}

template<typename TypeTag>
std::pair<std::optional<typename GasLiftStage2<TypeTag>::GradPairItr>,
          std::optional<typename GasLiftStage2<TypeTag>::GradPairItr>>
GasLiftStage2<TypeTag>::OptimizeState::
getEcoGradients(std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    if (inc_grads.size() > 0 && dec_grads.size() > 0) {
        sortGradients_(inc_grads);
        sortGradients_(dec_grads);
        auto inc_grad = std::prev(inc_grads.end());
        std::optional<GradPairItr> inc_grad_opt;
        std::optional<GradPairItr> dec_grad_opt;
        for (auto itr = dec_grads.begin(); itr != dec_grads.end(); itr++) {
            if (itr->first == inc_grad->first) {
                // Don't consider decremental gradients with the same well name
                continue;
            }
            dec_grad_opt = itr;
            break;
        }
        if (dec_grad_opt) {
            inc_grad_opt = inc_grad;
            return { inc_grad_opt, dec_grad_opt };
        }
    }
    return {std::nullopt, std::nullopt};
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
recalculateGradients(
         std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads,
         GradPairItr &min_dec_grad_itr, GradPairItr &max_inc_grad_itr)
{
    {
        const std::string &name = max_inc_grad_itr->first;
        GasLiftSingleWell &gs_well = *(this->parent.stage1_wells_.at(name).get());
        auto inc_grad = calcIncOrDecGrad_(name, gs_well, /*increase=*/true);
        if (inc_grad) {
            max_inc_grad_itr->second = inc_grad->grad;
            auto old_grad = this->parent.updateIncGrad_(name, *inc_grad);
            // The old incremental gradient becomes the new decremental gradient
            this->parent.updateDecGrad_(name, old_grad);
        }
        else {
            inc_grads.erase(max_inc_grad_itr); // NOTE: this invalidates max_inc_grad_itr
            auto old_grad = this->parent.deleteIncGradItem_(name);
            // NOTE: Either creates a new item or reassigns
            this->parent.dec_grads_[name] = old_grad;
        }
    }
    {
        const std::string &name = min_dec_grad_itr->first;
        GasLiftSingleWell &gs_well = *(this->parent.stage1_wells_.at(name).get());
        auto dec_grad = calcIncOrDecGrad_(name, gs_well, /*increase=*/false);
        if (dec_grad) {
            min_dec_grad_itr->second = dec_grad->grad;
            auto old_grad = this->parent.updateDecGrad_(name, *dec_grad);
            // The old decremental gradient becomes the new incremental gradient
            this->parent.updateIncGrad_(name, old_grad);
        }
        else {
            dec_grads.erase(min_dec_grad_itr); // NOTE: this invalidates min_dec_grad_itr
            auto old_grad = this->parent.deleteDecGradItem_(name);
            // NOTE: Either creates a new item or reassigns
            this->parent.inc_grads_[name] = old_grad;
        }
    }
}

// Take one ALQ increment from well1, and give it to well2
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
    redistributeALQ( GradPairItr &min_dec_grad, GradPairItr &max_inc_grad)
{
    const std::string msg = fmt::format(
        "redistributing ALQ from well {} (dec gradient: {}) "
        "to well {} (inc gradient {})",
        min_dec_grad->first, min_dec_grad->second,
        max_inc_grad->first, max_inc_grad->second);
    displayDebugMessage_(msg);
    addOrRemoveALQincrement_(
        this->parent.dec_grads_, /*well_name=*/min_dec_grad->first, /*add=*/false);
    addOrRemoveALQincrement_(
        this->parent.inc_grads_, /*well_name=*/max_inc_grad->first, /*add=*/true);
}

/****************************************
 * Private methods
 ****************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
addOrRemoveALQincrement_(GradMap &grad_map, const std::string well_name, bool add)
{
    GasLiftSingleWell &gs_well = *(this->parent.stage1_wells_.at(well_name).get());
    const GradInfo &gi = grad_map.at(well_name);
    gs_well.addOrSubtractAlqIncrement(gi, add);
}

template<typename TypeTag>
std::optional<typename GasLiftStage2<TypeTag>::GradInfo>
GasLiftStage2<TypeTag>::OptimizeState::
calcIncOrDecGrad_(
    const std::string well_name, const GasLiftSingleWell &gs_well, bool increase)
{
    /*
    {
        const std::string msg = fmt::format("well: {}", well_name);
        displayDebugMessage_(msg);
    }
    */
    auto grad = gs_well.calcIncOrDecGradient(increase);
    if (grad) {
        const std::string msg = fmt::format(
            "well {} : adding {} gradient = {}",
            well_name,
            (increase ? "incremental" : "decremental"),
            grad->grad
        );
        displayDebugMessage_(msg);
    }
    else {
        /*
        const std::string msg = fmt::format(
            "well {} : not able to obtain {} gradient",
            well_name,
            (increase ? "incremental" : "decremental")
        );
        displayDebugMessage_(msg);
        */
    }
    return grad;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
displayDebugMessage_(const std::string &msg)
{
    this->parent.displayDebugMessage_(msg, this->group.name());
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
displayWarning_(const std::string &msg)
{
    this->parent.displayWarning_(msg, this->group.name());
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
sortGradients_(std::vector<GradPair> &grads)
{
    auto cmp = [](GradPair a, GradPair b) {
         return a.second <  b.second;
    };
    std::sort(grads.begin(), grads.end(), cmp);
}

