// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Test for efficiency factor handling in checkGroupConstraintsProd()
 */

#define BOOST_TEST_MODULE EfficiencyFactor

#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>

#include <boost/test/unit_test.hpp>

namespace Opm::Properties::TTag {
    struct TestEfficiencyFactorTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

//! \brief Test fixture with unit conversion helpers for metric units
struct MetricUnitFixture
{
    //! Convert SI rate (m³/s) to metric rate (SM3/day)
    static double metric_rate(double si_rate)
    {
        using namespace Opm::unit;
        return convert::to(si_rate, cubic(meter) / day);
    }
};

} // anonymous namespace

BOOST_FIXTURE_TEST_SUITE(EfficiencyFactorTests, MetricUnitFixture)

BOOST_AUTO_TEST_CASE(TestEfficiencyFactorUnits)
{
    using TypeTag = Opm::Properties::TTag::TestEfficiencyFactorTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;

    const std::string filename = "EFFICIENCY_FACTOR.DATA";

    auto simulator = Opm::initSimulator<TypeTag>(filename.data());

    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(86400);  // 1 day in seconds
    simulator->model().newtonMethod().setIterationIndex(0);

    WellModel& well_model = simulator->problem().wellModel();
    const auto& schedule = simulator->vanguard().schedule();
    const int report_step_idx = 0;

    // Initialize the well model
    // - beginReportStep(): sets up report step data
    // - beginTimeStep(): calls updateAndCommunicateGroupData() which computes
    //   reduction rates via updateGroupTargetReduction()
    //
    // Note: We intentionally do NOT call updateWellControls() because it
    // internally calls checkGroupConstraintsProd(), which is the method we want
    // to test in isolation.
    well_model.beginReportStep(report_step_idx);
    well_model.beginTimeStep();

    Opm::DeferredLogger deferred_logger;

    // Get the group hierarchy info
    const auto& plat_group = schedule.getGroup("PLAT", report_step_idx);
    const auto& mani_group = schedule.getGroup("MANI", report_step_idx);

    // Verify efficiency factors are set correctly
    const double E_PLAT = plat_group.getGroupEfficiencyFactor();
    const double E_MANI = mani_group.getGroupEfficiencyFactor();

    BOOST_CHECK_CLOSE(E_PLAT, 0.9, 1e-8);
    BOOST_CHECK_CLOSE(E_MANI, 0.8, 1e-8);

    // Get current rates from group state
    const auto& group_state = well_model.groupState();
    auto& group_state_helper = well_model.groupStateHelper();
    const auto& pu = well_model.phaseUsage();

    const int water_phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
    const int gas_phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
    const int oil_phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);

    // First, let's check the group hierarchy as read from Schedule
    std::cout << "\n=== Group Hierarchy ===" << std::endl;
    for (const std::string group_name : {"FIELD", "PLAT", "MANI", "MANI2"}) {
        const auto& grp = schedule.getGroup(group_name, report_step_idx);
        std::cout << group_name << ":" << std::endl;
        std::cout << "  child groups: ";
        for (const auto& child : grp.groups()) {
            std::cout << child << " ";
        }
        std::cout << (grp.groups().empty() ? "(none)" : "") << std::endl;
        std::cout << "  wells: ";
        for (const auto& well : grp.wells()) {
            std::cout << well << " ";
        }
        std::cout << (grp.wells().empty() ? "(none)" : "") << std::endl;
    }

    // production_rates() returns total rates for a group including all descendants
    // This is different from sumWellSurfaceRates() which only includes direct wells
    std::cout << "\n=== Group Production Rates (from GroupState) ===" << std::endl;
    for (const std::string group_name : {"FIELD", "PLAT", "MANI", "MANI2"}) {
        if (group_state.has_production_rates(group_name)) {
            const auto& rates = group_state.production_rates(group_name);
            std::cout << "production_rates(" << group_name << ", oil): "
                      << metric_rate(rates[oil_phase_pos]) << " SM3/day" << std::endl;
        } else {
            std::cout << "production_rates(" << group_name << "): not set" << std::endl;
        }
    }


    // Get sumWellSurfaceRates for MANI
    // This is what's used to compute current_rate_available
    // Use proper phase index lookup instead of hardcoding
    const double sum_rates_mani = group_state_helper.sumWellSurfaceRates(
        mani_group, oil_phase_pos, /*is_injector=*/false);
    BOOST_CHECK_CLOSE(metric_rate(sum_rates_mani), 5000.0, 1e-2);

    std::cout << "sumWellSurfaceRates(MANI, oil): " << metric_rate(sum_rates_mani)
              << " SM3/day" << std::endl;

    // Get MANI2 info
    const auto& mani2_group = schedule.getGroup("MANI2", report_step_idx);
    const double E_MANI2 = mani2_group.getGroupEfficiencyFactor();
    std::cout << "MANI2 efficiency (E_MANI2): " << E_MANI2 << std::endl;

    const double sum_rates_mani2 = group_state_helper.sumWellSurfaceRates(
        mani2_group, oil_phase_pos, /*is_injector=*/false);
    std::cout << "sumWellSurfaceRates(MANI2, oil): " << metric_rate(sum_rates_mani2)
              << " SM3/day" << std::endl;

    // Get individual well rates from well state
    const auto& well_state = well_model.wellState();
    std::cout << "\n=== Individual Well Rates ===" << std::endl;
    for (const std::string well_name : {"WELL-A", "WELL-B"}) {
        const auto well_idx = well_state.index(well_name);
        if (well_idx.has_value()) {
            const auto& ws = well_state.well(well_idx.value());
            std::cout << well_name << " oil rate: "
                      << metric_rate(ws.surface_rates[oil_phase_pos]) << " SM3/day" << std::endl;
            std::cout << "  production_cmode: " << static_cast<int>(ws.production_cmode)
                      << " (GRUP=" << static_cast<int>(Opm::Well::ProducerCMode::GRUP) << ")" << std::endl;
        }
    }

    // Check groupControlledWells for both groups
    const int mani_num_gr_ctrl = group_state_helper.groupControlledWells(
        "MANI", "", /*is_producer=*/true, Opm::Phase::OIL);
    std::cout << "groupControlledWells(MANI): " << mani_num_gr_ctrl << std::endl;

    // Check if PLAT has any direct wells
    const auto& plat_wells = plat_group.wells();
    std::cout << "PLAT direct wells: " << plat_wells.size() << std::endl;
    for (const auto& well_name : plat_wells) {
        std::cout << "  " << well_name << std::endl;
    }

    // Get reduction rates for all groups
    std::cout << "\n=== Group Reduction Rates ===" << std::endl;
    for (const std::string group_name : {"FIELD", "PLAT", "MANI", "MANI2"}) {
        if (group_state.has_production_reduction_rates(group_name)) {
            const auto& reduction = group_state.production_reduction_rates(group_name);
            std::cout << group_name << " reduction (oil): "
                      << metric_rate(reduction[oil_phase_pos]) << " SM3/day" << std::endl;
        } else {
            std::cout << group_name << " reduction: not set" << std::endl;
        }
    }

    // Manual verification of PLAT reduction
    // According to updateGroupTargetReductionRecursive_():
    // - MANI has individual control → PLAT_reduction += E_MANI * sumWellSurfaceRates(MANI)
    // - MANI2 has FLD control, no guide rate → PLAT_reduction += E_MANI2 * MANI2_sub_reduction
    std::cout << "\n=== Manual Verification ===" << std::endl;
    std::cout << "MANI contribution to PLAT reduction:" << std::endl;
    std::cout << "  = E_MANI * sumWellSurfaceRates(MANI)" << std::endl;
    std::cout << "  = " << E_MANI << " * " << metric_rate(sum_rates_mani) << std::endl;
    std::cout << "  = " << metric_rate(E_MANI * sum_rates_mani) << " SM3/day" << std::endl;

    // For MANI2: since it has FLD control and no guide rate, its sub_group_target_reduction
    // is accumulated. But sub_group_target_reduction is the value RETURNED from the recursive
    // call, not the stored production_reduction_rates.
    //
    // Looking at updateGroupTargetReductionRecursive_():
    // - For MANI2 (FLD control, num_group_controlled_wells > 0, no guide rate):
    //   → Uses sub_group_target_reduction (returned from recursion)
    // - But what does MANI2's recursion return?
    //   → WELL-B is under GRUP control, so it's NOT added to MANI2's reduction
    //   → MANI2's sub_group_target_reduction should be 0
    //
    // Wait - let's check WELL-B's control mode. If WELL-B is under GRUP control,
    // it doesn't add to the reduction. But sumWellSurfaceRates still includes it!
    //
    // Actually, looking at lines 1771-1779:
    //   if (ws.production_cmode != Well::ProducerCMode::GRUP)
    //       group_target_reduction[phase] -= ws.surface_rates[phase] * efficiency;
    //
    // So wells under GRUP control do NOT contribute to reduction!
    // But wait - MANI2 has FLD control with num_group_controlled_wells > 0,
    // so it takes the "else" branch at line 1727, NOT the sumWellSurfaceRates path.
    //
    // The sub_group_target_reduction for MANI2 should include WELL-B only if
    // WELL-B is NOT under GRUP control... but WELL-B IS under GRUP control.
    //
    // So why is PLAT reduction = 7201? Let's check if there's something else going on.

    std::cout << "MANI2 contribution to PLAT reduction:" << std::endl;
    std::cout << "  MANI2 has FLD control, no guide rate" << std::endl;
    std::cout << "  → Code path: uses sub_group_target_reduction (returned from recursion)" << std::endl;
    std::cout << "  → But WELL-B is under GRUP control, so it shouldn't be in MANI2's reduction" << std::endl;
    std::cout << std::endl;

    // Let's check what sumWellSurfaceRates returns vs the stored reduction
    std::cout << "  sumWellSurfaceRates(MANI2): " << metric_rate(sum_rates_mani2)
              << " SM3/day" << std::endl;

    double mani2_reduction_oil = 0.0;
    if (group_state.has_production_reduction_rates("MANI2")) {
        const auto& mani2_reduction = group_state.production_reduction_rates("MANI2");
        mani2_reduction_oil = mani2_reduction[oil_phase_pos];
        std::cout << "  production_reduction_rates(MANI2): "
                  << metric_rate(mani2_reduction_oil) << " SM3/day" << std::endl;
    }

    // The actual MANI2 contribution uses sumWellSurfaceRates because:
    // - MANI2 has FLD control (individual_control = false)
    // - num_group_controlled_wells for MANI2 > 0 (WELL-B under GRUP)
    // - MANI2 has no guide rate
    // So code takes the path at lines 1727-1736: sub_group_efficiency * sub_group_target_reduction
    //
    // But sub_group_target_reduction for MANI2 comes from the recursion INTO MANI2,
    // which only adds wells NOT under GRUP control. Since WELL-B IS under GRUP,
    // it shouldn't be added... unless there's something we're missing.
    //
    // Actually wait - let me re-read the code. The recursion computes MANI2's
    // sub_group_target_reduction, but WELL-B under GRUP means it's NOT added
    // to that reduction at line 1775. So MANI2's sub_group_target_reduction = 0.
    //
    // But then where does the 7201 come from?
    //
    // Oh! Maybe the issue is that MANI2's num_group_controlled_wells = 0?
    // If MANI2 has no group-controlled wells, it would take the path at line 1722,
    // which uses sumWellSurfaceRates instead!

    // Check MANI2's num_group_controlled_wells
    const int mani2_num_gr_ctrl = group_state_helper.groupControlledWells(
        "MANI2", "", /*is_producer=*/true, Opm::Phase::OIL);
    std::cout << "  groupControlledWells(MANI2): " << mani2_num_gr_ctrl << std::endl;

    // Check PLAT's num_group_controlled_wells (needed for local_reduction_level)
    const int plat_num_gr_ctrl = group_state_helper.groupControlledWells(
        "PLAT", "", /*is_producer=*/true, Opm::Phase::OIL);
    std::cout << "groupControlledWells(PLAT): " << plat_num_gr_ctrl << std::endl;

    // If mani2_num_gr_ctrl == 0, then MANI2's contribution uses sumWellSurfaceRates
    double mani2_contribution = 0.0;
    if (mani2_num_gr_ctrl == 0) {
        // Path at line 1722: sub_group_efficiency * sumWellSurfaceRates
        mani2_contribution = E_MANI2 * sum_rates_mani2;
        std::cout << "  → Since num_group_controlled_wells = 0, uses sumWellSurfaceRates path" << std::endl;
        std::cout << "  MANI2 contribution = E_MANI2 * sumWellSurfaceRates(MANI2)" << std::endl;
        std::cout << "                     = " << E_MANI2 << " * " << metric_rate(sum_rates_mani2) << std::endl;
        std::cout << "                     = " << metric_rate(mani2_contribution) << " SM3/day" << std::endl;
    } else {
        // Path at line 1727: sub_group_efficiency * sub_group_target_reduction
        mani2_contribution = E_MANI2 * mani2_reduction_oil;
        std::cout << "  → Since num_group_controlled_wells > 0, uses sub_group_target_reduction path" << std::endl;
        std::cout << "  MANI2 contribution = E_MANI2 * MANI2_sub_group_target_reduction" << std::endl;
        std::cout << "                     = " << E_MANI2 << " * " << metric_rate(mani2_reduction_oil) << std::endl;
        std::cout << "                     = " << metric_rate(mani2_contribution) << " SM3/day" << std::endl;
    }

    // Total expected PLAT reduction
    double expected_plat_reduction = E_MANI * sum_rates_mani + mani2_contribution;
    std::cout << "\nExpected PLAT reduction (oil):" << std::endl;
    std::cout << "  = MANI_contribution + MANI2_contribution" << std::endl;
    std::cout << "  = " << metric_rate(E_MANI * sum_rates_mani) << " + "
              << metric_rate(mani2_contribution) << std::endl;
    std::cout << "  = " << metric_rate(expected_plat_reduction) << " SM3/day" << std::endl;

    // Additional debugging: Check if sumWellSurfaceRates(PLAT) includes both wells
    // sumWellSurfaceRates recursively sums rates with efficiency factors
    const double sum_rates_plat = group_state_helper.sumWellSurfaceRates(
        plat_group, oil_phase_pos, /*is_injector=*/false);
    std::cout << "\nsumWellSurfaceRates(PLAT, oil): " << metric_rate(sum_rates_plat)
              << " SM3/day" << std::endl;
    // Expected well rates in SM3/day: WELL-A ~5001, WELL-B ~4443
    const double well_a_rate = 5001.0 / 86400.0;  // ~0.0578704 m3/s
    const double well_b_rate = 4443.0 / 86400.0;  // ~0.0514194 m3/s
    std::cout << "  Expected: E_MANI * WELL-A + E_MANI2 * WELL-B" << std::endl;
    std::cout << "          = " << E_MANI << " * " << metric_rate(well_a_rate) << " + "
              << E_MANI2 << " * " << metric_rate(well_b_rate) << std::endl;
    double expected_sum_plat = E_MANI * well_a_rate + E_MANI2 * well_b_rate;
    std::cout << "          = " << metric_rate(expected_sum_plat) << " SM3/day" << std::endl;
    std::cout << "  Note: sumWellSurfaceRates includes efficiency factors of subgroups!" << std::endl;

    // Compare with actual PLAT reduction
    if (group_state.has_production_reduction_rates("PLAT")) {
        const auto& plat_reduction = group_state.production_reduction_rates("PLAT");
        double actual_plat_reduction = plat_reduction[oil_phase_pos];
        std::cout << "\nActual PLAT reduction (oil): " << metric_rate(actual_plat_reduction)
                  << " SM3/day" << std::endl;
        std::cout << "Difference: " << metric_rate(actual_plat_reduction - expected_plat_reduction)
                  << " SM3/day" << std::endl;
    }

    // The key test: verify that sumWellSurfaceRates includes/excludes
    // appropriate efficiency factors
    //
    // According to our analysis:
    // - sumWellSurfaceRates(MANI) includes E_child factors (none for single well)
    //   but NOT E_MANI itself
    // - This becomes current_rate_available in checkGroupConstraintsProd()
    //
    // The question is: when we divide target by efficiency_factor (E_PLAT * E_MANI),
    // does target_rate_available end up in the same units as current_rate_available?

    std::cout << "\n=== Group Control Status ===" << std::endl;
    if (group_state.has_production_control("MANI")) {
        auto ctrl = group_state.production_control("MANI");
        std::cout << "MANI control mode: " << Opm::Group::ProductionCMode2String(ctrl) << std::endl;
        std::cout << "  individual_control(MANI) = "
                  << ((ctrl != Opm::Group::ProductionCMode::FLD && ctrl != Opm::Group::ProductionCMode::NONE) ? "true" : "false")
                  << std::endl;
    }
    if (group_state.has_production_control("MANI2")) {
        auto ctrl = group_state.production_control("MANI2");
        std::cout << "MANI2 control mode: " << Opm::Group::ProductionCMode2String(ctrl) << std::endl;
        std::cout << "  individual_control(MANI2) = "
                  << ((ctrl != Opm::Group::ProductionCMode::FLD && ctrl != Opm::Group::ProductionCMode::NONE) ? "true" : "false")
                  << std::endl;
    }
    if (group_state.has_production_control("PLAT")) {
        auto ctrl = group_state.production_control("PLAT");
        std::cout << "PLAT control mode: " << Opm::Group::ProductionCMode2String(ctrl) << std::endl;
    }

    // Verify guide rate is set for PLAT and MANI2
    const auto& guide_rate = well_model.guideRate();
    bool plat_has_guide_rate = guide_rate.has("PLAT");
    std::cout << "PLAT has guide rate: " << (plat_has_guide_rate ? "yes" : "no") << std::endl;
    bool mani2_has_guide_rate = guide_rate.has("MANI2");
    std::cout << "MANI2 has guide rate: " << (mani2_has_guide_rate ? "yes" : "no") << std::endl;

    // Bug 002 verification: Call checkGroupConstraintsProd directly with known rates
    //
    // Test scenario (updated):
    // - FIELD has ORAT control with target 10000 sm3/day
    // - PLAT has FLD control with guide rate (FORM)
    // - MANI has individual ORAT control with target 5000 sm3/day (triggers check)
    // - MANI2 has FLD control (makes groupControlledWells(PLAT) > 0)
    // - efficiency_factor = E_MANI * E_PLAT = 0.8 * 0.9 = 0.72
    //
    // When calling checkGroupConstraintsProd for MANI checking against PLAT:
    // - rates = sumWellSurfaceRates(MANI) (includes WEFAC, excludes E_MANI)
    // - current_rate_available = -calcModeRateFromRates(rates)
    // - target_rate_available = target / efficiency_factor
    //
    // Since PLAT has FLD control with guide rate, it will recurse to FIELD.
    // local_reduction_level should be set at PLAT level (where guide rate is)
    // because groupControlledWells(PLAT) > 0 (due to MANI2/WELL-B under FLD/GRUP).
    //
    // Bug 002 question: Are current_rate_available and target_rate_available in the same units?

    std::cout << "\n=== Bug 002 Verification ===" << std::endl;

    // Create a mock rate vector representing MANI's surface rates
    // These would come from sumWellSurfaceRates in the real call
    // For testing, we use a known value
    const int num_phases = pu.numActivePhases();
    std::vector<double> mock_rates(num_phases, 0.0);

    // Set oil rate to 3000 sm3/day (negative for production)
    // Convert to SI units (m³/s) since the simulator works internally in SI
    // 1 day = 86400 seconds
    const double rate_sm3_per_day = 3000.0;
    const double seconds_per_day = 86400.0;
    mock_rates[oil_phase_pos] = -rate_sm3_per_day / seconds_per_day;  // ≈ -0.0347 m³/s

    // Setup resv_coeff (conversion factors, just use 1.0 for this test)
    std::vector<double> resv_coeff(num_phases, 1.0);

    // Initial efficiency factor starts as E_MANI (0.8)
    // Will be multiplied by E_PLAT (0.9) when recursing through PLAT
    const double initial_efficiency = E_MANI;

    std::cout << "Test scenario:" << std::endl;
    std::cout << "  FIELD target (ORAT): 10000 SM3/day" << std::endl;
    std::cout << "  PLAT control: FLD (with guide rate)" << std::endl;
    std::cout << "  MANI control: ORAT (5000 SM3/day) - triggers checkGroupConstraintsProd" << std::endl;
    std::cout << "  MANI2 control: FLD - makes groupControlledWells(PLAT) > 0" << std::endl;
    std::cout << "  Mock MANI oil rate: " << rate_sm3_per_day << " SM3/day" << std::endl;
    std::cout << "  E_MANI: " << E_MANI << std::endl;
    std::cout << "  E_PLAT: " << E_PLAT << std::endl;
    std::cout << "  Expected accumulated efficiency: " << (E_MANI * E_PLAT) << std::endl;
    std::cout << std::endl;

    // Make the actual call to checkGroupConstraintsProd
    // This is the function we want to debug to understand Bug 002
    //
    // Call signature:
    // checkGroupConstraintsProd(name, parent, group, rates, efficiency_factor, resv_coeff, check_guide_rate, deferred_logger)
    //
    // We call it for MANI checking against its parent PLAT
    // Since PLAT has FLD control with guide rate, it will recurse to FIELD
    // The efficiency_factor starts as E_MANI (0.8) and will be multiplied by E_PLAT (0.9)

    std::cout << "Calling checkGroupConstraintsProd:" << std::endl;
    std::cout << "  name: MANI" << std::endl;
    std::cout << "  parent: PLAT" << std::endl;
    std::cout << "  group: PLAT (will recurse to FIELD since PLAT has FLD control)" << std::endl;
    std::cout << "  rates[oil]: " << metric_rate(mock_rates[oil_phase_pos]) << " SM3/day" << std::endl;
    std::cout << "  initial efficiency_factor: " << initial_efficiency << " (E_MANI)" << std::endl;
    std::cout << "  (will become " << (E_MANI * E_PLAT) << " after recursing through PLAT)" << std::endl;
    std::cout << std::endl;

    // NOTE: Set a breakpoint here in gdb to step through checkGroupConstraintsProd
    // Key variables to watch:
    // - efficiency_factor (should accumulate to 0.72)
    // - local_reduction_level (should be > 0 at PLAT level where guide rate is)
    // - current_rate_available (from rates parameter)
    // - target_rate_available (target / efficiency_factor at line 413)
    auto [constraint_violated, scale] = group_state_helper.checkGroupConstraintsProd(
        "MANI",                    // name (the group being checked)
        "PLAT",                    // parent (immediate parent group)
        plat_group,                // group (PLAT - will recurse to FIELD)
        mock_rates.data(),         // rates
        initial_efficiency,        // efficiency_factor (E_MANI, will accumulate E_PLAT)
        resv_coeff,                // resv_coeff
        /*check_guide_rate=*/true, // check_guide_rate
        deferred_logger            // deferred_logger
    );

    std::cout << "checkGroupConstraintsProd result:" << std::endl;
    std::cout << "  constraint_violated: " << (constraint_violated ? "yes" : "no") << std::endl;
    std::cout << "  scale: " << scale << std::endl;
    std::cout << std::endl;

    // Bug 002 Analysis:
    // The key question is whether current_rate_available and target_rate_available
    // are in the same units when compared at line 419
    std::cout << "Bug 002 Analysis:" << std::endl;
    std::cout << "  current_rate_available source: rates parameter (mock_rates)" << std::endl;
    std::cout << "  - In real usage: sumWellSurfaceRates(MANI) - reduction_rates" << std::endl;
    std::cout << "  - Includes: WEFAC (well efficiency)" << std::endl;
    std::cout << "  - Excludes: E_MANI (MANI's own efficiency)" << std::endl;
    std::cout << std::endl;
    std::cout << "  target_rate_available = target / efficiency_factor" << std::endl;
    std::cout << "  - efficiency_factor should be " << (E_MANI * E_PLAT) << std::endl;
    std::cout << "  - This converts target to 'full rate' units" << std::endl;
    std::cout << std::endl;
    std::cout << "  Use gdb to step through and verify the intermediate values." << std::endl;
    std::cout << "  Watch for local_reduction_level > 0 when PLAT has guide rate." << std::endl;

    // Verify the key setup conditions for Bug 002 scenario
    BOOST_CHECK(plat_has_guide_rate);  // PLAT needs guide rate for local_reduction_level
    BOOST_CHECK_CLOSE(E_PLAT * E_MANI, 0.72, 1e-6);  // Verify accumulated efficiency

    well_model.endTimeStep();
    well_model.endEpisode();

    std::cout << "\n=== Test Complete ===" << std::endl;
    std::cout << "Test setup verified for Bug 002 scenario." << std::endl;
    std::cout << "Manual inspection or debug logging needed to verify intermediate values." << std::endl;
}

BOOST_AUTO_TEST_CASE(TestIsolatedEfficiencyFactor)
{
    //
    // This test creates an isolated environment with manually set well rates
    // and group states, avoiding the complexity of beginTimeStep() which uses
    // nupcolWellState and various update functions that obscure the test setup.
    //
    using TypeTag = Opm::Properties::TTag::TestEfficiencyFactorTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;

    const std::string filename = "EFFICIENCY_FACTOR.DATA";
    auto simulator = Opm::initSimulator<TypeTag>(filename.data());

    // Basic simulator setup (no beginReportStep/beginTimeStep)
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(86400);  // 1 day in seconds
    simulator->model().newtonMethod().setIterationIndex(0);

    WellModel& well_model = simulator->problem().wellModel();
    const auto& schedule = simulator->vanguard().schedule();
    const int report_step_idx = 0;
    const auto& pu = well_model.phaseUsage();

    const int oil_phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);

    // ========================================================================
    // Initialize wells via beginReportStep
    // ========================================================================
    // We need beginReportStep() to create the wells and initialize WellState.
    // This sets up wells from the Schedule but doesn't call updateAndCommunicateGroupData().
    well_model.beginReportStep(report_step_idx);

    // ========================================================================
    // STEP 1: Set up fixed well rates in WellState
    // ========================================================================
    // We want predictable rates:
    //   WELL-A: 5000 SM3/day oil (negative for production)
    //   WELL-B: 4000 SM3/day oil (negative for production)

    const double well_a_oil_rate = -5000.0 / 86400.0;  // m³/s (negative = production)
    const double well_b_oil_rate = -4000.0 / 86400.0;  // m³/s

    auto& ws = well_model.wellState();

    // Set WELL-A rates
    {
        auto well_idx = ws.index("WELL-A");
        BOOST_REQUIRE(well_idx.has_value());
        auto& well_a = ws.well(well_idx.value());
        std::fill(well_a.surface_rates.begin(), well_a.surface_rates.end(), 0.0);
        well_a.surface_rates[oil_phase_pos] = well_a_oil_rate;
        well_a.production_cmode = Opm::Well::ProducerCMode::GRUP;  // Under group control
    }

    // Set WELL-B rates
    {
        auto well_idx = ws.index("WELL-B");
        BOOST_REQUIRE(well_idx.has_value());
        auto& well_b = ws.well(well_idx.value());
        std::fill(well_b.surface_rates.begin(), well_b.surface_rates.end(), 0.0);
        well_b.surface_rates[oil_phase_pos] = well_b_oil_rate;
        well_b.production_cmode = Opm::Well::ProducerCMode::GRUP;  // Under group control
    }

    // ========================================================================
    // STEP 2: nupcolWellState
    // ========================================================================
    // Note: updateNupcolWGState() is protected, but we don't need it since
    // sumWellSurfaceRates uses the active wellState(), not nupcolWellState.

    // ========================================================================
    // STEP 3: Set up group control modes in GroupState
    // ========================================================================
    auto& gs = well_model.groupState();

    // MANI has individual ORAT control (this triggers checkGroupConstraintsProd)
    gs.production_control("MANI", Opm::Group::ProductionCMode::ORAT);

    // MANI2 has FLD control (follows parent)
    gs.production_control("MANI2", Opm::Group::ProductionCMode::FLD);

    // PLAT has FLD control (follows parent FIELD)
    gs.production_control("PLAT", Opm::Group::ProductionCMode::FLD);

    // FIELD has ORAT control
    gs.production_control("FIELD", Opm::Group::ProductionCMode::ORAT);

    // ========================================================================
    // STEP 4: Update group-controlled wells count
    // ========================================================================
    Opm::DeferredLogger deferred_logger;
    auto& gsh = well_model.groupStateHelper();

    // Update the groupControlledWells count based on well control modes
    gsh.updateGroupControlledWells(/*is_production_group=*/true, Opm::Phase::OIL, deferred_logger);

    // ========================================================================
    // STEP 5: Update group target reduction rates
    // ========================================================================
    const auto& field_group = schedule.getGroup("FIELD", report_step_idx);
    gsh.updateGroupTargetReduction(field_group, /*is_injector=*/false);

    // ========================================================================
    // VERIFY: Check the setup
    // ========================================================================
    std::cout << "\n=== Isolated Test Setup ===" << std::endl;

    // Verify well rates
    std::cout << "\nWell rates (from WellState):" << std::endl;
    std::cout << "  WELL-A oil: " << metric_rate(ws.well("WELL-A").surface_rates[oil_phase_pos])
              << " SM3/day" << std::endl;
    std::cout << "  WELL-B oil: " << metric_rate(ws.well("WELL-B").surface_rates[oil_phase_pos])
              << " SM3/day" << std::endl;

    // Verify sumWellSurfaceRates matches our manual setup
    const auto& mani_group = schedule.getGroup("MANI", report_step_idx);
    const auto& mani2_group = schedule.getGroup("MANI2", report_step_idx);
    const auto& plat_group = schedule.getGroup("PLAT", report_step_idx);

    const double sum_mani = gsh.sumWellSurfaceRates(mani_group, oil_phase_pos, /*injector=*/false);
    const double sum_mani2 = gsh.sumWellSurfaceRates(mani2_group, oil_phase_pos, /*injector=*/false);
    const double sum_plat = gsh.sumWellSurfaceRates(plat_group, oil_phase_pos, /*injector=*/false);

    std::cout << "\nsumWellSurfaceRates:" << std::endl;
    std::cout << "  MANI:  " << metric_rate(sum_mani) << " SM3/day (expected: 5000)" << std::endl;
    std::cout << "  MANI2: " << metric_rate(sum_mani2) << " SM3/day (expected: 4000)" << std::endl;
    std::cout << "  PLAT:  " << metric_rate(sum_plat) << " SM3/day" << std::endl;

    // Verify efficiency factors
    const double E_MANI = mani_group.getGroupEfficiencyFactor();
    const double E_MANI2 = mani2_group.getGroupEfficiencyFactor();
    const double E_PLAT = plat_group.getGroupEfficiencyFactor();

    std::cout << "\nEfficiency factors:" << std::endl;
    std::cout << "  E_MANI:  " << E_MANI << std::endl;
    std::cout << "  E_MANI2: " << E_MANI2 << std::endl;
    std::cout << "  E_PLAT:  " << E_PLAT << std::endl;

    // Expected PLAT sumWellSurfaceRates with efficiency factors
    const double expected_sum_plat = E_MANI * (-well_a_oil_rate) + E_MANI2 * (-well_b_oil_rate);
    std::cout << "\nExpected sumWellSurfaceRates(PLAT):" << std::endl;
    std::cout << "  = E_MANI * WELL-A + E_MANI2 * WELL-B" << std::endl;
    std::cout << "  = " << E_MANI << " * 5000 + " << E_MANI2 << " * 4000" << std::endl;
    std::cout << "  = " << metric_rate(expected_sum_plat) << " SM3/day" << std::endl;

    // Verify group reduction rates
    std::cout << "\nGroup reduction rates:" << std::endl;
    for (const std::string gname : {"FIELD", "PLAT", "MANI", "MANI2"}) {
        if (gs.has_production_reduction_rates(gname)) {
            const auto& red = gs.production_reduction_rates(gname);
            std::cout << "  " << gname << ": " << metric_rate(red[oil_phase_pos]) << " SM3/day" << std::endl;
        } else {
            std::cout << "  " << gname << ": not set" << std::endl;
        }
    }

    // Verify groupControlledWells
    std::cout << "\ngroupControlledWells:" << std::endl;
    std::cout << "  MANI:  " << gsh.groupControlledWells("MANI", "", true, Opm::Phase::OIL) << std::endl;
    std::cout << "  MANI2: " << gsh.groupControlledWells("MANI2", "", true, Opm::Phase::OIL) << std::endl;
    std::cout << "  PLAT:  " << gsh.groupControlledWells("PLAT", "", true, Opm::Phase::OIL) << std::endl;

    // ========================================================================
    // TEST: Verify sumWellSurfaceRates returns expected values
    // ========================================================================
    BOOST_CHECK_CLOSE(metric_rate(sum_mani), 5000.0, 1e-2);
    BOOST_CHECK_CLOSE(metric_rate(sum_mani2), 4000.0, 1e-2);

    // PLAT should sum child groups with their efficiency factors
    // sumWellSurfaceRates(PLAT) = E_MANI * WELL-A + E_MANI2 * WELL-B
    //                          = 0.8 * 5000 + 0.7 * 4000 = 4000 + 2800 = 6800
    BOOST_CHECK_CLOSE(metric_rate(sum_plat), 6800.0, 1e-2);

    std::cout << "\n=== Isolated Test Complete ===" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
