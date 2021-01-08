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
    WellState &well_state,
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
    schedule_{ebos_simulator.vanguard().schedule()},
    phase_usage_{well_model_.phaseUsage()},
    glo_{schedule_.glo(report_step_idx_)},
    debug_{true}
{ }

/********************************************
 * Public methods in alphabetical order
 ********************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
runOptimize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);

    optimizeGroupsRecursive_(group);

}


/********************************************
 * Private methods in alphabetical order
 ********************************************/

// Update "stage2_state" in GasLiftSingleWell for "well_name" to the
//   new ALQ value and related data (the data has already been computed and
//   saved in "grad_map")
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
addOrRemoveALQincrement_(GradMap &grad_map, const std::string well_name, bool add)
{
    GasLiftSingleWell &gs_well = *(this->stage1_wells_.at(well_name).get());
    const GradInfo &gi = grad_map.at(well_name);
    gs_well.updateStage2State(gi, add);
    this->well_state_.setALQ(well_name, gi.alq);
}

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
displayWarning_(const std::string &msg)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION (STAGE2) : {}", msg);
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

template<typename TypeTag>
std::tuple<double, double, double>
GasLiftStage2<TypeTag>::
getCurrentGroupRates_(const Opm::Group &group)
{
    double oil_rate = 0.0;
    double gas_rate = 0.0;
    double alq = 0.0;
    // NOTE: A group can either contain wells or groups, but not both
    if (group.wellgroup()) {
        for (const std::string& well_name : group.wells()) {
            auto [sw_oil_rate, sw_gas_rate, sw_alq] =
                getCurrentWellRates_(well_name, group.name());
            oil_rate += sw_oil_rate;
            gas_rate += sw_gas_rate;
            alq += sw_alq;
        }
    }
    else {
        for (const std::string& group_name : group.groups()) {
            if(this->schedule_.hasGroup(group_name)) {
                const Group& sub_group =
                    this->schedule_.getGroup(group_name, this->report_step_idx_);
                // If groups have efficiency factors to model
                // synchronized downtime of their subordinate wells
                // (see keyword GEFAC), their lift gas injection rates
                // are multiplied by their efficiency factors when
                // they are added to the lift gas supply rate of the
                // parent group.
                const auto gefac = sub_group.getGroupEfficiencyFactor();
                auto [sg_oil_rate, sg_gas_rate, sg_alq] =
                    getCurrentGroupRates_(sub_group);
                oil_rate += (gefac * sg_oil_rate);
                gas_rate += (gefac * sg_gas_rate);
                alq += (gefac * sg_alq);
            }
        }
    }
    return std::make_tuple(oil_rate, gas_rate, alq);
}

template<typename TypeTag>
std::tuple<double, double, double>
GasLiftStage2<TypeTag>::
getCurrentWellRates_(const std::string &well_name, const std::string &group_name)
{
    double oil_rate, gas_rate, alq;
    bool success = false;
    const WellInterface<TypeTag> *well_ptr = nullptr;
    std::string debug_info;
    if (this->stage1_wells_.count(well_name) == 1) {
        GasLiftSingleWell &gs_well = *(this->stage1_wells_.at(well_name).get());
        const WellInterface<TypeTag> &well = gs_well.getStdWell();
        well_ptr = &well;
        if (gs_well.hasStage2Rates()) {
            std::tie(oil_rate, gas_rate) = gs_well.getStage2Rates();
            success = true;
            if (this->debug_) debug_info = "(A)";
        }
        else {
            std::tie(oil_rate, gas_rate) = getStdWellRates_(well);
            success = true;
            if (this->debug_) debug_info = "(B)";
        }
    }
    if (this->prod_wells_.count(well_name) == 1) {
        well_ptr = this->prod_wells_.at(well_name);
        std::tie(oil_rate, gas_rate) = getStdWellRates_(*well_ptr);
        success = true;
        if ( this->debug_) debug_info = "(C)";
    }
    if (success) {
        assert(well_ptr);
        assert(well_ptr->isProducer());
        alq = this->well_state_.getALQ(well_name);
        if (this->debug_) {
            const std::string msg = fmt::format(
                "Rates {} for well {} : oil: {}, gas: {}, alq: {}",
                debug_info, well_name, oil_rate, gas_rate, alq);
            displayDebugMessage_(msg, group_name);
        }
        // If wells have efficiency factors to take account of regular
        // downtime (see keyword WEFAC), their lift gas injection
        // rates are multiplied by their efficiency factors when they
        // are added to the group lift gas supply rate. This is
        // consistent with the summation of flow rates for wells with
        // downtime, and preserves the ratio of production rate to
        // lift gas injection rate.
        const auto &well_ecl = well_ptr->wellEcl();
        double factor = well_ecl.getEfficiencyFactor();
        oil_rate *= factor;
        gas_rate *= factor;
        alq *= factor;
        if (this->debug_ && (factor != 1)) {
            const std::string msg = fmt::format(
                "Well {} : efficiency factor {}. New rates : oil: {}, gas: {}, alq: {}",
                well_name, factor, oil_rate, gas_rate, alq);
            displayDebugMessage_(msg, group_name);
        }
    }
    else {
        // NOTE: This happens for wells that are not producers, or not active.
        if (this->debug_) {
            const std::string msg = fmt::format("Could not determine current rates for "
                "well {}: (not active or injector)", well_name);
            displayDebugMessage_(msg, group_name);
        }
        oil_rate = 0.0; gas_rate = 0.0; alq = 0.0;
    }
    return std::make_tuple(oil_rate, gas_rate, alq);
}

template<typename TypeTag>
std::pair<double, double>
GasLiftStage2<TypeTag>::
getStdWellRates_(const WellInterface<TypeTag> &well)
{
    const int well_index = well.indexOfWell();
    const int np = this->well_state_.numPhases();
    const auto& pu = well.phaseUsage();
    auto oil_rate =
        -this->well_state_.wellRates()[np * well_index + pu.phase_pos[Oil]];
    auto gas_rate =
        -this->well_state_.wellRates()[np * well_index + pu.phase_pos[Gas]];
    return {oil_rate, gas_rate};
}

// Find all subordinate wells of a given group.
//
// NOTE: A group can either contain wells or groups, not both.
//   If it contains groups, we have to traverse those recursively to find the wells.
//
// NOTE: This means that wells are located at the leaf nodes of the tree, and
//       groups are located at the other nodes (not leaf nodes) of the tree
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
        if(this->schedule_.hasGroup(group_name)) {
            const Group& sub_group =
                this->schedule_.getGroup(group_name, this->report_step_idx_);
            getGroupGliftWellsRecursive_(sub_group, wells);
        }
    }
    for (const std::string& well_name : group.wells()) {
        if (this->stage1_wells_.count(well_name) == 1) {
            GasLiftSingleWell *well_ptr = this->stage1_wells_.at(well_name).get();
            wells.push_back(well_ptr);
        }
    }
}

// NOTE: This method is called by optimizeGroupsRecursive_() such that for the
//   example group tree:
//
//                                       FIELD
//                                         |
//                                       PLAT-A
//                          ---------------+---------------------
//                         |                                    |
//                        M5S                                  M5N
//                ---------+----------                     -----+-------
//               |                   |                    |            |
//              B1                  G1                   C1           F1
//           ----+------          ---+---              ---+---       ---+---
//          |    |     |         |      |             |      |      |      |
//        B-1H  B-2H  B-3H     G-3H    G-4H         C-1H   C-2H    F-1H   F-2H
//
//

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
optimizeGroup_(const Opm::Group &group)
{
    const auto &gl_group = this->glo_.group(group.name());
    const auto &max_glift = gl_group.max_lift_gas();
    const auto &max_total_gas = gl_group.max_total_gas();
    if (group.has_control(Group::ProductionCMode::ORAT)
                       || max_glift || max_total_gas)
    {
        displayDebugMessage_("optimizing", group.name());
        auto wells = getGroupGliftWells_(group);
        std::vector<GradPair> inc_grads;
        std::vector<GradPair> dec_grads;
        redistributeALQ_(wells, group, inc_grads, dec_grads);
        removeSurplusALQ_(wells, group, dec_grads);
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
        if(!this->schedule_.hasGroup(group_name))
            continue;
        const Group& sub_group = this->schedule_.getGroup(
            group_name, this->report_step_idx_);
        optimizeGroupsRecursive_(sub_group);
    }
    // TODO: should we also optimize groups that do not have GLIFTOPT defined?
    //   (i.e. glo_.has_group(name) returns true)
    //   IF GLIFTOPT is not specified for the group or if item 2 of GLIFTOPT
    //   is defaulted, there is no maximum lift gas supply for the group.
    //   But even if there is no limit on the liftgas supply it can still
    //   be desireable to use as little ALQ as possible to achieve a
    //   group oil rate limit or gas rate limit.
    if (this->glo_.has_group(group.name())) // only optimize if GLIFTOPT is given
        optimizeGroup_(group);

}

// TODO: maybe the redistribution can be simplified by only doing it for the
//   FIELD group: For example:
//
//
//
//                                       FIELD
//                                         |
//                                       PLAT-A
//                          ---------------+---------------------
//                         |                                    |
//                        M5S                                  M5N
//                ---------+----------                     -----+-------
//               |                   |                    |            |
//              B1                  G1                   C1           F1
//           ----+------          ---+---              ---+---       ---+---
//          |    |     |         |      |             |      |      |      |
//        B-1H  B-2H  B-3H     G-3H    G-4H         C-1H   C-2H    F-1H   F-2H
//
//
//  it is probably unecessary to first redistribute ALQ for the wells B-1H, B-2H,
//  and B-3H in group B1, and then in a next separate step, redistribute ALQ for the
//  wells G-3H, and G-4H, and similarly, for the wells in groups C1, and F1,
//  and then, for the wells under M5S, and then M5N, and finally repeat the procedure
//  for all the wells under group PLAT-A, i.e. the wells:
//    B-1H, B-2H, B-3H, G-3H, G-4H, C-1H, C-2H, F-1H, and F-2H.
//  It seems like it would be more efficient to
//  just do it once for the topmost group "PLAT-A" and then skip redistribution for
//  all sub groups of "PLAT-A"
//
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
redistributeALQ_(std::vector<GasLiftSingleWell *> &wells,  const Opm::Group &group,
    std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    OptimizeState state {*this, group};
    if (!state.checkAtLeastTwoWells(wells)) {
        return;
    }
    assert(wells.size() >= 2);
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
            // Redistribute if the largest incremental gradient exceeds the
            //   smallest decremental gradient
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

// The group has surplus lift gas if it exceeds any production rate limits
//   or a lift gas supply limit, or contains any wells that have a weighted
//   decremental gradient less than the minimum economic gradient.
// Lift gas increments are removed in turn from the well that currently has
//   the smallest weighted decremental gradient, until there is no surplus
//   lift gas in the group.
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
removeSurplusALQ_(std::vector<GasLiftSingleWell *> &wells,  const Opm::Group &group,
    std::vector<GradPair> &dec_grads)
{
    if (dec_grads.size() == 0) {
        displayDebugMessage_("remove surplus ALQ: no wells to remove from");
        return;
    }
    assert(dec_grads.size() > 0);
    const auto &gl_group = this->glo_.group(group.name());
    const auto &max_glift = gl_group.max_lift_gas();
    const auto controls = group.productionControls(this->summary_state_);
    //const auto &max_total_gas = gl_group.max_total_gas();
    auto [oil_rate, gas_rate, alq] = getCurrentGroupRates_(group);
    auto min_eco_grad = this->glo_.min_eco_gradient();
    bool stop_iteration = false;
    while (!stop_iteration) {
            if (dec_grads.size() >= 2) {
                sortGradients_(dec_grads);
            }
            bool remove = false;  // remove ALQ in this iteration?
            auto dec_grad_itr = dec_grads.begin();
            if (dec_grad_itr->second < min_eco_grad) {
                remove = true;
            }
            else if (group.has_control(Group::ProductionCMode::ORAT)) {
                if (controls.oil_target < oil_rate  ) {
                    remove = true;
                }
            }
            else if (group.has_control(Group::ProductionCMode::GRAT)) {
                if (controls.gas_target < gas_rate  ) {
                    remove = true;
                }
            }
            else if (max_glift) {
                if (*max_glift < alq) {
                    remove = true;
                }
            }
            if (remove) {
                const auto &name = dec_grad_itr->first;
                addOrRemoveALQincrement_( this->dec_grads_, name, /*add=*/false);
                GasLiftSingleWell &gs_well = *(this->stage1_wells_.at(name).get());
                auto grad = gs_well.calcIncOrDecGradient(/*increase=*/false);
                // NOTE: We only update the decremental gradient, we
                //   do not need to update the incremental gradient,
                //   since it will no longer be used for anything..
                if (grad) {
                    dec_grad_itr->second = grad->grad;
                    updateDecGrad_(name, *grad);
                }
                else {
                    // NOTE: this invalidates min_dec_grad_itr
                    dec_grads.erase(dec_grad_itr);
                    auto old_grad = this->parent.deleteDecGradItem_(name);
                    if (dec_grads.size() == 0) stop_iteration = true;
                }
            }
            else {
                stop_iteration = true;
            }
        }
    }
    //std::vector<GradPair> dec_grads;
    //dec_grads.reserve(wells.size());
    //state.calculateEcoGradients(wells, inc_grads, dec_grads);
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
void
GasLiftStage2<TypeTag>::
sortGradients_(std::vector<GradPair> &grads)
{
    auto cmp = [](GradPair a, GradPair b) {
         return a.second <  b.second;
    };
    std::sort(grads.begin(), grads.end(), cmp);
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

/***********************************************
 * Public methods declared in OptimizeState
 ***********************************************/

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
        this->parent.sortGradients_(inc_grads);
        this->parent.sortGradients_(dec_grads);
        // The largest incremental gradient is the last element
        auto inc_grad = std::prev(inc_grads.end());
        std::optional<GradPairItr> inc_grad_opt;
        std::optional<GradPairItr> dec_grad_opt;
        // The smallest decremental gradient is at the beginning
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
            return { dec_grad_opt, inc_grad_opt };
        }
    }
    return {std::nullopt, std::nullopt};
}

// Recalculate gradients (and related information, see struct GradInfo in
//   GasLiftSingleWell.hpp) after an ALQ increment
//   has been given from the well with minumum decremental gradient (represented
//   by the input argument min_dec_grad_itr) to the well with the largest
//   incremental gradient (represented by input argument max_inc_grad_itr).
//
// For the well with the largest incremental gradient, we compute a new
//   incremental gradient given the new ALQ. The new decremental gradient for this
//   well is set equal to the current incremental gradient (before the ALQ is added)
// Similiarly, for the well with the smalles decremental gradient, we compute
//   a new decremental gradient given the new ALQ. The new incremental gradient
//   for this well is set equal to the current decremental gradient
//   (before the ALQ is subtracted)
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
            updateGradVector_(name, dec_grads, old_grad.grad);
        }
        else {
            inc_grads.erase(max_inc_grad_itr); // NOTE: this invalidates max_inc_grad_itr
            auto old_grad = this->parent.deleteIncGradItem_(name);
            // NOTE: Either creates a new item or reassigns
            this->parent.dec_grads_[name] = old_grad;
            updateGradVector_(name, dec_grads, old_grad.grad);
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
            updateGradVector_(name, inc_grads, old_grad.grad);
        }
        else {
            dec_grads.erase(min_dec_grad_itr); // NOTE: this invalidates min_dec_grad_itr
            auto old_grad = this->parent.deleteDecGradItem_(name);
            // NOTE: Either creates a new item or reassigns
            this->parent.inc_grads_[name] = old_grad;
            updateGradVector_(name, inc_grads, old_grad.grad);
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
    this->parent.addOrRemoveALQincrement_(
        this->parent.dec_grads_, /*well_name=*/min_dec_grad->first, /*add=*/false);
    this->parent.addOrRemoveALQincrement_(
        this->parent.inc_grads_, /*well_name=*/max_inc_grad->first, /*add=*/true);
}

/**********************************************
 * Private methods declared in OptimizeState
 **********************************************/

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
updateGradVector_(const std::string &name, std::vector<GradPair> &grads, double grad)
{
    for (auto itr = grads.begin(); itr != grads.end(); itr++) {
        if (itr->first == name) {
            itr->second = grad;
            return;
        }
    }
    grads.push_back({name, grad});
    // NOTE: the gradient vector is no longer sorted, but sorting will be done
    //   later in getEcoGradients()
}
