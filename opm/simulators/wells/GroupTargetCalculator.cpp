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
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/simulators/wells/GroupTargetCalculator.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <tuple>  // for std::tie
#include <fmt/format.h>

namespace Opm {

// ----------------------------------------------------
// Constructor for the GroupTargetCalculator class
// ----------------------------------------------------
template<class Scalar, class IndexTraits>
GroupTargetCalculator<Scalar, IndexTraits>::
GroupTargetCalculator(
    const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
    const GroupStateHelperType& group_state_helper,
    DeferredLogger& deferred_logger
) :
    well_model_{well_model},
    group_state_helper_{group_state_helper},
    well_state_{group_state_helper.wellState()},
    group_state_{group_state_helper.groupState()},
    schedule_{group_state_helper.schedule()},
    summary_state_{group_state_helper.summaryState()},
    phase_usage_{group_state_helper.phaseUsage()},
    guide_rate_{group_state_helper.guideRate()},
    report_step_idx_{group_state_helper.reportStepIdx()},
    deferred_logger_{deferred_logger},
    resv_coeffs_inj_(group_state_helper.phaseUsage().numPhases, 0.0)
{
    std::tie(this->fipnum_, this->pvtreg_) = this->well_model_.getGroupFipnumAndPvtreg();
    // Get the reservoir coefficients for injection
    // NOTE: The production reservoir coefficients depend on the group production rates, so we
    //  cannot precompute them here.
    this->well_model_.calcInjResvCoeff(this->fipnum_, this->pvtreg_, this->resv_coeffs_inj_);
}

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::InjectionTargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
groupInjectionTarget(const Group& group, ReservoirCoupling::Phase injection_phase)
{
    GeneralCalculator general_calculator{*this, group, injection_phase};
    auto target_info = general_calculator.calculateGroupTarget();
    if (!target_info) {
        // No target was calculated, return an empty optional.
        return std::nullopt;
    }
    return std::make_optional<InjectionTargetInfo>(
        InjectionTargetInfo{target_info->target, std::get<Group::InjectionCMode>(target_info->cmode)}
    );
}

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::ProductionTargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
groupProductionTarget(const Group& group) {
    GeneralCalculator general_calculator{*this, group};
    auto target_info = general_calculator.calculateGroupTarget();
    if (!target_info) {
        // No target was calculated, return an empty optional.
        return std::nullopt;
    }
    return std::make_optional<ProductionTargetInfo>(
        ProductionTargetInfo{
            target_info->target,
            std::get<Group::ProductionCMode>(target_info->cmode)
        }
    );
}

// ----------------------------------------------------
// Constructor for inner class GeneralCalculator
// ----------------------------------------------------

template<class Scalar, class IndexTraits>
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
GeneralCalculator(
    GroupTargetCalculator<Scalar, IndexTraits>& parent_calculator,
    const Group& original_group,
    std::optional<ReservoirCoupling::Phase> injection_phase  // Only used for injectors
) :
    parent_calculator_{parent_calculator},
    original_group_{original_group},
    injection_phase_{injection_phase},
    resv_coeffs_prod_(this->phaseUsage().numPhases, 0.0)
{
    this->wellModel().calcResvCoeff(
        this->fipnum(),
        this->pvtreg(),
        this->groupState().production_rates(original_group.name()),
        this->resv_coeffs_prod_
    );
}

// -------------------------------------------------------
// Public methods for the GeneralCalculator class
// -------------------------------------------------------

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
calculateGroupTarget()
{
    // TODO: For now this is adapted to the reservoir coupling implementation where "group" below will
    //   correspond to a master group. In the future we might want to port the logic to include
    //   e.g. checkGroupConstraintsProd() and checkGroupConstraintsInj() in GroupStateHelper.cpp.
    const auto& group = this->original_group_;  // The bottom group we want to calculate the target for.
    if (group.is_field() || !this->parentGroupControlAvailable_(group)) {
        return this->getTargetNoGuideRate_(group);
    }
    assert(this->parentGroupControlAvailable_(group));
    if (!this->hasGuideRate_(group)) {
        if (this->hasFldOrNoneControl_(group)) {
            // Parent control is available, but no guide rate is defined. This is illegal for a master group
            // under FLD or NONE control.
            throw std::runtime_error("No guide rate defined for master group " + group.name());
        }
        // A master group with:
        //   - not FLD or NONE control,
        //   - parent control available,
        //   - but no guide rate
        // this could be considered an error, but we can also fall back to use its own target.
        return this->getTargetNoGuideRate_(group);
    }
    const auto efficiency_factor = group.getGroupEfficiencyFactor();
    auto target = this->calculateGroupTargetRecursive_(group, efficiency_factor);
    if (target) {
        // TODO: Should we switch group to FLD control mode now?
        return target;
    }
    // If a higher level target was not found or not violated, use the group's own target.
    return this->getTargetNoGuideRate_(group);
}

template<class Scalar, class IndexTraits>
void
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
defLogThrow(const std::string& message)
{
    OPM_DEFLOG_THROW(std::runtime_error, message, this->deferredLogger());
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getGratSalesInjectionTarget(const Group& group) const
{
    const auto& gconsale = this->gconsale();
    const auto& group_name = group.name();
    if (gconsale.has(group_name)) {
        return gconsale.get(group_name, this->summaryState()).sales_target;
    }
    return 0.0;

}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getGratSalesProductionTarget(const Group& group) const
{
    if (this->groupState().has_grat_sales_target(group.name()))
        return this->groupState().grat_sales_target(group.name());
    return 0.0;
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::GeneralCalculator::TargetCalculatorType
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getInjectionTargetCalculator(const Group& group)
{
    const auto& control_mode = this->groupState().injection_control(group.name(), this->injectionPhase_());
    Scalar sales_target = this->getGratSalesInjectionTarget(group);
    bool has_gpmaint_control = group.has_gpmaint_control(this->injectionPhase_(), control_mode);
    return InjectionTargetCalculator{
        control_mode,
        this->phaseUsage(),
        this->resvCoeffsInj(),
        group.name(),
        sales_target,
        this->groupState(),
        this->injectionPhase_(),
        has_gpmaint_control,
        this->deferredLogger()
    };
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::GeneralCalculator::TargetCalculatorType
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getProductionTargetCalculator(const Group& group) const
{
    const auto& control_mode = this->groupState().production_control(group.name());
    Scalar grat_sales_target = this->getGratSalesProductionTarget(group);
    bool has_gpmaint_control = group.has_gpmaint_control(control_mode);
    return TargetCalculator{
        control_mode,
        this->phaseUsage(),
        this->resvCoeffsProd(),
        grat_sales_target,
        group.name(),
        this->groupState(),
        has_gpmaint_control
    };
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::GeneralCalculator::TargetCalculatorType
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getTargetCalculator(const Group& group)
{
    if (this->isInjector()) {
        return this->getInjectionTargetCalculator(group);
    }
    else {
        return this->getProductionTargetCalculator(group);
    }
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getTargetFromCalculator(const TargetCalculatorType& target_calculator, const Group& group)
{
    if (this->isInjector()) {
        const auto& control_mode = this->groupState().injection_control(group.name(), this->injectionPhase_());
        std::optional<Group::InjectionControls> ctrl;
        if (!group.has_gpmaint_control(this->injectionPhase_(), control_mode))
            ctrl = group.injectionControls(this->injectionPhase_(), this->summaryState());
        return TargetInfo{
            std::get<InjectionTargetCalculator>(target_calculator).groupTarget(
                ctrl, this->deferredLogger()
            ), control_mode
        };
    }
    else {
        const auto& control_mode = this->groupState().production_control(group.name());
        std::optional<Group::ProductionControls> ctrl;
        if (!group.has_gpmaint_control(control_mode))
            ctrl = group.productionControls(this->summaryState());
        return TargetInfo{
            std::get<TargetCalculator>(target_calculator).groupTarget(
                ctrl, this->deferredLogger()
            ), control_mode
        };
    }
}


// -------------------------------------------------------
// Private methods for the GeneralCalculator class
// -------------------------------------------------------

// See comments above RescoupTargetCalculator::calculateMasterGroupTargetsAndSendToSlaves() about how
// the target is calculated.
template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
calculateGroupTargetRecursive_(const Group& group, const Scalar efficiency_factor)
{
    if (this->hasFldOrNoneControl_(group)) {
        // NOTE: A pure injection group is assumed to have production control NONE and
        //   a pure production group is assumed to have injection control NONE, see
        //   GroupStateHelper.cpp::setCmodeGroup() for details.
        if (!this->parentGroupControlAvailable_(group)) {
            // If we are here, the group has control mode NONE, and either
            //   - item 8 of GCONPROD/GCONINJE is "NO" or,
            //   - group.is_field() is true.
            // NOTE: it cannot have control mode FLD, since FLD will ovverride item 8 of
            //   GCONPROD/GCONINJE, see GroupKeywordHandlers.cpp
            // TODO: We could throw an error here. Or try to fall back to use the group's own target.
            // Fall back to use the group's own target.
            return std::nullopt;
        } else {
            assert(!group.is_field());  // The field group should not have parent control.
            auto new_efficiency_factor = efficiency_factor * group.getGroupEfficiencyFactor();
            return this->calculateGroupTargetRecursive_(this->parentGroup(group), new_efficiency_factor);
        }
    }

    if ((this->isInjector() && !group.isInjectionGroup()) ||
            (this->isProducer() && !group.isProductionGroup()))
    {
        // Controlling group that is not of the type we are interested in.
        // This should never happen, since any group that is not an injector should have injection
        // control NONE. And any group that is not a producer should have production control NONE.
        // See GroupStateHelper.cpp::setCmodeGroup() for details.
        // .. and therefore should be caught by the hasFldOrNoneControl_() check above.
        this->defLogThrow(
            fmt::format("Controlling group that is not of the type we are interested in: {}", group.name())
        );
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    const auto& top_group = group;
    const auto& bottom_group = this->original_group_;
    TopToBottomCalculator top_bottom_calc{*this, top_group, bottom_group, efficiency_factor};
    return top_bottom_calc.calculateGroupTarget();
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getTargetNoGuideRate_(const Group& group)
{
    const auto target_calculator = this->getTargetCalculator(group);
    return this->getTargetFromCalculator(target_calculator, group);
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
hasFldOrNoneControl_(const Group& group)
{
    // Check if the group (constrained to production or injection) has control FLD or NONE.
    // For example, a pure injection group will have production control NONE.
    const auto& name = group.name();
    if (this->isInjector()) {
        return this->groupState().has_field_or_none_control(name, this->injectionPhase_());
    }
    else {
        return this->groupState().has_field_or_none_control(name);
    }
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
hasGuideRate_(const Group& group) const
{
    return this->guideRate().has(group.name());
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
parentGroupControlAvailable_(const Group& group)
{
    if (this->isInjector()) {
        return group.injectionGroupControlAvailable(this->injectionPhase_());
    } else {
        return group.productionGroupControlAvailable();
    }
}

template<class Scalar, class IndexTraits>
Phase
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
injectionPhase_()
{
    switch (this->injection_phase_.value()) {
        case ReservoirCoupling::Phase::Water:
            return Phase::WATER;
        case ReservoirCoupling::Phase::Oil:
            return Phase::OIL;
        case ReservoirCoupling::Phase::Gas:
            return Phase::GAS;
        default:
           this->defLogThrow("Invalid injection phase");
           return Phase::WATER; // Unreachable, but satisfies compiler
    }
}

// ----------------------------------------------------
// Constructor for inner class TopToBottomCalculator
// ----------------------------------------------------

template<class Scalar, class IndexTraits>
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
TopToBottomCalculator(
    GeneralCalculator& parent_calculator,
    const Group& top_group,
    const Group& bottom_group,
    Scalar efficiency_factor
) :
    parent_calculator_{parent_calculator},
    top_group_{top_group},
    bottom_group_{bottom_group},
    efficiency_factor_{efficiency_factor}
{
    if (this->isInjector()) {
        this->initForInjector_();
    }
    else {
        this->initForProducer_();
    }
}

// -------------------------------------------------------
// Public methods for the TopToBottomCalculator class
// -------------------------------------------------------

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
calculateGroupTarget()
{
    auto orig_target = this->getTopLevelTarget_();
    const auto chain = this->getGroupChainTopBot_();
    // Because the bottom group (the original group) is the last of the elements,
    //   and not an ancestor, we subtract one:
    const std::size_t num_ancestors = chain.size() - 1;
    Scalar target = orig_target;
    for (std::size_t i = 0; i < num_ancestors; ++i) {
        if ((i == 0) || this->guideRate().has(chain[i])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= this->localReduction_(chain[i]);
        }
        target *= this->localFraction_(chain[i + 1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    target = std::max(Scalar(0.0), target / this->efficiency_factor_);
    return TargetInfo{target, this->toplevel_control_mode_};
}

// -------------------------------------------------------
// Private methods for the TopToBottomCalculator class
// -------------------------------------------------------

template<class Scalar, class IndexTraits>
void
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
initForInjector_()
{
    auto calc = this->getInjectionTargetCalculator_(this->top_group_);
    auto guide_target_mode =
        std::get<InjectionTargetCalculator>(this->target_calculator_).guideTargetMode();
    this->fraction_calculator_.emplace(
        this->schedule(),
        this->groupStateHelper(),
        this->summaryState(),
        this->reportStepIdx(),
        &this->guideRate(),
        GuideRateModel::Target{guide_target_mode},
        /*is_producer=*/false,
        this->injectionPhase_()
    );
}

template<class Scalar, class IndexTraits>
void
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
initForProducer_()
{
    auto calc = this->getProductionTargetCalculator_(this->top_group_);
    this->target_calculator_.template emplace<TargetCalculator>(std::get<TargetCalculator>(calc));
    auto guide_target_mode =
        std::get<TargetCalculator>(this->target_calculator_).guideTargetMode();
    auto dummy_phase = Phase::OIL; // Dummy phase, not used for producers.
    this->fraction_calculator_.emplace(
        this->schedule(),
        this->groupStateHelper(),
        this->summaryState(),
        this->reportStepIdx(),
        &this->guideRate(),
        GuideRateModel::Target{guide_target_mode},
        /*is_producer=*/true,
        dummy_phase
    );
}

template<class Scalar, class IndexTraits>
std::vector<std::string>
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getGroupChainTopBot_() const
{
    return this->groupStateHelper().groupChainTopBot(this->bottom_group_.name(), this->top_group_.name());
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getTopLevelTarget_()
{
    if (this->isInjector()) {
        auto control_mode = this->groupState().injection_control(
            this->top_group_.name(), this->injectionPhase_()
        );
        std::optional<Group::InjectionControls> ctrl;
        if (!this->top_group_.has_gpmaint_control(this->injectionPhase_(), control_mode))
            ctrl = this->top_group_.injectionControls(this->injectionPhase_(), this->summaryState());
        return std::get<InjectionTargetCalculator>(this->target_calculator_).groupTarget(
            ctrl, this->deferredLogger()
        );
    }
    else {
        auto control_mode = this->groupState().production_control(this->top_group_.name());
        std::optional<Group::ProductionControls> ctrl;
        if (!this->top_group_.has_gpmaint_control(control_mode))
            ctrl = this->top_group_.productionControls(this->summaryState());
        return std::get<TargetCalculator>(this->target_calculator_).groupTarget(
            ctrl, this->deferredLogger()
        );
    }
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
localFraction_(const std::string& group_name)
{
    const auto& always_included_child = this->bottom_group_.name();
    return this->fraction_calculator_->localFraction(group_name, always_included_child);
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
localReduction_(const std::string& group_name)
{
    if (this->isInjector()) {
        const std::vector<Scalar>& group_target_reductions =
            this->groupState().injection_reduction_rates(group_name);
        return std::get<InjectionTargetCalculator>(this->target_calculator_).calcModeRateFromRates(
            group_target_reductions
        );
    }
    else {
        const std::vector<Scalar>& group_target_reductions =
            this->groupState().production_reduction_rates(group_name);
        return std::get<TargetCalculator>(this->target_calculator_).calcModeRateFromRates(
            group_target_reductions
        );

    }
}


template class GroupTargetCalculator<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class GroupTargetCalculator<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
