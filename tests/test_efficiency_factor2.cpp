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
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
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
    // STEP 2: Update global well info (isProductionGrup)
    // ========================================================================
    // After setting well control modes, we need to update the global well info
    // so that isProductionGrup() returns true for wells under group control.
    // This is needed for FractionCalculator::guideRateSum() to include wells.
    std::vector<Opm::Well::Status> well_status(ws.size(), Opm::Well::Status::OPEN);
    const auto& comm = simulator->vanguard().grid().comm();
    ws.updateGlobalIsGrup(comm, well_status);

    // ========================================================================
    // STEP 3: Mock guide rates
    // ========================================================================
    // Set up guide rates for wells to control how rates are distributed.
    // Without this, localFraction() falls back to potentials or returns 1.0.
    // We use guide rates proportional to well rates for simplicity.
    //
    // Guide rate values are in output units (SM3/day for metric), not SI units.
    // The init_grvalue() method stores values directly without conversion.
    auto& guide_rate = well_model.guideRate();
    const double sim_time = 0.0;

    // Set guide rates for wells (values in SM3/day, same as their production rates)
    guide_rate.init_grvalue(report_step_idx, "WELL-A",
        Opm::GuideRate::GuideRateValue{sim_time, 5000.0, Opm::GuideRateModel::Target::OIL});
    guide_rate.init_grvalue(report_step_idx, "WELL-B",
        Opm::GuideRate::GuideRateValue{sim_time, 4000.0, Opm::GuideRateModel::Target::OIL});

    // Set guide rate for PLAT group - this is needed to trigger local_reduction_level
    // at PLAT level in checkGroupConstraintsProd(). Without this, the algorithm
    // doesn't set local_reduction_level and the addback logic at line 407 is skipped.
    // The value here represents PLAT's total guide rate capacity.
    guide_rate.init_grvalue(report_step_idx, "PLAT",
        Opm::GuideRate::GuideRateValue{sim_time, 9000.0, Opm::GuideRateModel::Target::OIL});

    // ========================================================================
    // STEP 4: Set up group production rates in GroupState
    // ========================================================================
    // The checkGroupConstraintsProd function uses production_rates() to compute
    // current_rate_fraction when checking guide rate violations. We need to
    // populate these for the groups in the chain.
    auto& gs = well_model.groupState();
    const int num_phases = pu.numActivePhases();

    // Get group objects and their efficiency factors
    const auto& mani_group = schedule.getGroup("MANI", report_step_idx);
    const auto& mani2_group = schedule.getGroup("MANI2", report_step_idx);
    const auto& plat_group = schedule.getGroup("PLAT", report_step_idx);
    const auto& field_group = schedule.getGroup("FIELD", report_step_idx);

    const double E_MANI = mani_group.getGroupEfficiencyFactor();
    const double E_MANI2 = mani2_group.getGroupEfficiencyFactor();
    const double E_PLAT = plat_group.getGroupEfficiencyFactor();

    // Set production rates for groups (these are the "current" group rates)
    // MANI produces 5000 SM3/day, but with efficiency factor this shows up as:
    // MANI: 5000 SM3/day (direct well rate)
    // PLAT: E_MANI * 5000 + E_MANI2 * 4000 = 6800 SM3/day
    // FIELD: E_PLAT * 6800 = 0.9 * 6800 = 6120 SM3/day
    std::vector<double> mani_prod_rates(num_phases, 0.0);
    mani_prod_rates[oil_phase_pos] = -well_a_oil_rate;  // 5000 SM3/day (positive)

    std::vector<double> mani2_prod_rates(num_phases, 0.0);
    mani2_prod_rates[oil_phase_pos] = -well_b_oil_rate;  // 4000 SM3/day (positive)

    std::vector<double> plat_prod_rates(num_phases, 0.0);
    plat_prod_rates[oil_phase_pos] = E_MANI * (-well_a_oil_rate) + E_MANI2 * (-well_b_oil_rate);  // 6800 SM3/day

    std::vector<double> field_prod_rates(num_phases, 0.0);
    field_prod_rates[oil_phase_pos] = E_PLAT * plat_prod_rates[oil_phase_pos];  // 6120 SM3/day

    gs.update_production_rates("MANI", mani_prod_rates);
    gs.update_production_rates("MANI2", mani2_prod_rates);
    gs.update_production_rates("PLAT", plat_prod_rates);
    gs.update_production_rates("FIELD", field_prod_rates);

    // ========================================================================
    // STEP 5: Set up group control modes in GroupState
    // ========================================================================

    // MANI has individual ORAT control (this triggers checkGroupConstraintsProd)
    gs.production_control("MANI", Opm::Group::ProductionCMode::ORAT);

    // MANI2 has FLD control (follows parent)
    gs.production_control("MANI2", Opm::Group::ProductionCMode::FLD);

    // PLAT has FLD control (follows parent FIELD)
    gs.production_control("PLAT", Opm::Group::ProductionCMode::FLD);

    // FIELD has ORAT control
    gs.production_control("FIELD", Opm::Group::ProductionCMode::ORAT);

    // ========================================================================
    // STEP 5: Update group-controlled wells count
    // ========================================================================
    Opm::DeferredLogger deferred_logger;
    auto& gsh = well_model.groupStateHelper();

    // Update the groupControlledWells count based on well control modes
    gsh.updateGroupControlledWells(/*is_production_group=*/true, Opm::Phase::OIL, deferred_logger);

    // ========================================================================
    // STEP 6: Update group target reduction rates
    // ========================================================================
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
    const double sum_mani = gsh.sumWellSurfaceRates(mani_group, oil_phase_pos, /*injector=*/false);
    const double sum_mani2 = gsh.sumWellSurfaceRates(mani2_group, oil_phase_pos, /*injector=*/false);
    const double sum_plat = gsh.sumWellSurfaceRates(plat_group, oil_phase_pos, /*injector=*/false);

    std::cout << "\nsumWellSurfaceRates:" << std::endl;
    std::cout << "  MANI:  " << metric_rate(sum_mani) << " SM3/day (expected: 5000)" << std::endl;
    std::cout << "  MANI2: " << metric_rate(sum_mani2) << " SM3/day (expected: 4000)" << std::endl;
    std::cout << "  PLAT:  " << metric_rate(sum_plat) << " SM3/day" << std::endl;

    // Verify efficiency factors (already defined in STEP 4)
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

    // Verify isProductionGrup is set correctly
    std::cout << "\nisProductionGrup:" << std::endl;
    std::cout << "  WELL-A: " << std::boolalpha << ws.isProductionGrup("WELL-A") << std::endl;
    std::cout << "  WELL-B: " << std::boolalpha << ws.isProductionGrup("WELL-B") << std::endl;

    // Verify guide rates are set
    std::cout << "\nGuide rates:" << std::endl;
    std::cout << "  WELL-A has guide rate: " << std::boolalpha << guide_rate.has("WELL-A") << std::endl;
    std::cout << "  WELL-B has guide rate: " << std::boolalpha << guide_rate.has("WELL-B") << std::endl;
    std::cout << "  PLAT has guide rate: " << std::boolalpha << guide_rate.has("PLAT") << std::endl;

    // ========================================================================
    // TEST: Verify sumWellSurfaceRates returns expected values
    // ========================================================================
    BOOST_CHECK_CLOSE(metric_rate(sum_mani), 5000.0, 1e-2);
    BOOST_CHECK_CLOSE(metric_rate(sum_mani2), 4000.0, 1e-2);

    // PLAT should sum child groups with their efficiency factors
    // sumWellSurfaceRates(PLAT) = E_MANI * WELL-A + E_MANI2 * WELL-B
    //                          = 0.8 * 5000 + 0.7 * 4000 = 4000 + 2800 = 6800
    BOOST_CHECK_CLOSE(metric_rate(sum_plat), 6800.0, 1e-2);

    // ========================================================================
    // TEST: Call checkGroupConstraintsProd() for MANI
    // ========================================================================
    // This is the function we're investigating for Bug 002.
    //
    // Setup:
    // - MANI has individual ORAT control (5000 SM3/day limit in DATA file)
    // - MANI's parent is PLAT (FLD control)
    // - PLAT's parent is FIELD (ORAT control = 10000 SM3/day)
    // - Efficiency factors: E_MANI=0.8, E_PLAT=0.9
    //
    // The function will:
    // 1. See PLAT has FLD control → recurse to FIELD
    // 2. At FIELD, check ORAT constraint against MANI's rates
    // 3. Use accumulated efficiency_factor = E_MANI * E_PLAT = 0.72
    //
    std::cout << "\n=== checkGroupConstraintsProd Test ===" << std::endl;

    // Prepare rates array for MANI (negative for production)
    // This represents the "current rate" that will be compared against the target
    std::vector<double> mani_rates(num_phases, 0.0);
    mani_rates[oil_phase_pos] = well_a_oil_rate;  // -5000 SM3/day in SI units

    // Prepare resv_coeff (reservoir volume coefficients) - use 1.0 for simplicity
    std::vector<double> resv_coeff(num_phases, 1.0);

    // Initial efficiency factor is E_MANI (will accumulate E_PLAT during recursion)
    const double initial_efficiency = E_MANI;

    std::cout << "\nCalling checkGroupConstraintsProd for MANI:" << std::endl;
    std::cout << "  name: MANI (the entity being checked)" << std::endl;
    std::cout << "  parent: PLAT (immediate parent)" << std::endl;
    std::cout << "  group: PLAT (starting group for recursion)" << std::endl;
    std::cout << "  rates[oil]: " << metric_rate(mani_rates[oil_phase_pos]) << " SM3/day" << std::endl;
    std::cout << "  initial efficiency_factor: " << initial_efficiency << " (E_MANI)" << std::endl;
    std::cout << "  (will become " << (E_MANI * E_PLAT) << " after recursing through PLAT)" << std::endl;

    // Call checkGroupConstraintsProd
    // Arguments: name, parent, group, rates, efficiency_factor, resv_coeff, check_guide_rate, deferred_logger
    auto [constraint_violated, scale] = gsh.checkGroupConstraintsProd(
        "MANI",                    // name (the group/well being checked)
        "PLAT",                    // parent (immediate parent group name)
        plat_group,                // group (PLAT - will recurse to FIELD since PLAT has FLD)
        mani_rates.data(),         // rates (MANI's production rates)
        initial_efficiency,        // efficiency_factor (E_MANI, will accumulate E_PLAT)
        resv_coeff,                // resv_coeff
        /*check_guide_rate=*/true, // check_guide_rate
        deferred_logger            // deferred_logger
    );

    std::cout << "\nResult:" << std::endl;
    std::cout << "  constraint_violated: " << std::boolalpha << constraint_violated << std::endl;
    std::cout << "  scale: " << scale << std::endl;

    // Compute what the target_rate_available would have been
    // scale = target_rate_available / current_rate_available
    // current_rate_available = |mani_rates[oil]| = 5000 SM3/day
    const double current_rate_available = -mani_rates[oil_phase_pos];  // Make positive
    const double target_rate_available = scale * current_rate_available;
    std::cout << "\nDerived values:" << std::endl;
    std::cout << "  current_rate_available: " << metric_rate(current_rate_available) << " SM3/day" << std::endl;
    std::cout << "  target_rate_available: " << metric_rate(target_rate_available) << " SM3/day" << std::endl;
    std::cout << "  (target_rate_available = scale × current_rate_available)" << std::endl;

    // ========================================================================
    // Analysis of expected result
    // ========================================================================
    // FIELD has ORAT = 10000 SM3/day
    // MANI is producing 5000 SM3/day (WELL-A)
    // MANI2 is producing 4000 SM3/day (WELL-B) under FLD control
    //
    // At FIELD level, the target calculation involves:
    // - Starting target: 10000 SM3/day (FIELD's ORAT)
    // - Subtract PLAT's reduction (depends on MANI2's contribution)
    // - Potentially add back MANI's contribution (if local_reduction_level applies)
    // - Divide by efficiency_factor to get target_rate_available
    //
    // The constraint is violated if current_rate > target_rate_available
    //
    // With MANI producing 5000 SM3/day and FIELD allowing 10000 SM3/day total,
    // there should be headroom, so constraint_violated should be false.
    //
    std::cout << "\nExpected behavior:" << std::endl;
    std::cout << "  FIELD ORAT target: 10000 SM3/day" << std::endl;
    std::cout << "  MANI current rate: 5000 SM3/day" << std::endl;
    std::cout << "  If no constraint violation, scale > 1.0 (headroom available)" << std::endl;

    // ========================================================================
    // Manual calculation to verify the result
    // ========================================================================
    // The checkGroupConstraintsProd function computes:
    //   target_rate_available = target / efficiency_factor
    // where target is computed by iterating through the chain and applying
    // reductions and guide rate fractions.
    //
    // Chain: FIELD -> PLAT -> MANI
    // - At FIELD (ii=0): target starts at FIELD's ORAT = 10000 SM3/day
    // - At PLAT (ii=1): target -= local_reduction(PLAT), then *= local_fraction
    //
    // Let's compute what we expect:
    //
    // FIELD reduction = E_PLAT * PLAT_reduction = 0.9 * 4000 = 3600 SM3/day
    // (This matches our output above!)
    //
    // At FIELD level (ii=0):
    //   target = 10000 SM3/day (FIELD's ORAT)
    //   local_reduction = FIELD_reduction = 3600 SM3/day
    //   Since local_reduction_level >= 0, we DON'T subtract at this level
    //   target stays 10000
    //   local_fraction = 1.0 (assuming MANI gets all of FIELD's allocation)
    //   target = 10000 * 1.0 = 10000
    //
    // At PLAT level (ii=1):
    //   local_reduction = PLAT_reduction = 4000 SM3/day
    //   target -= 4000 -> target = 6000 SM3/day
    //
    //   If local_reduction_level == 1:
    //     Addback: target += current_rate_available * efficiency_factor
    //            = target + 5000 * 0.72 = target + 3600
    //            = 6000 + 3600 = 9600 SM3/day
    //
    //   local_fraction = guide_rate(MANI) / total_guide_rate
    //                  = 5000 / (5000 + 4000) = 5000/9000 = 0.556
    //   target *= 0.556 -> target = 9600 * 0.556 = 5333 SM3/day
    //
    // Final: target_rate_available = target / efficiency_factor
    //                              = 5333 / 0.72 = 7407 SM3/day
    //
    // But we got 8169.93 SM3/day, so the calculation differs.
    // This suggests the local_fraction or local_reduction values differ.

    std::cout << "\n=== Manual Calculation Verification ===" << std::endl;

    // Get FIELD's ORAT target
    const auto& field_group_ctrl = field_group.productionControls(simulator->vanguard().summaryState());
    const double field_orat_target = field_group_ctrl.oil_target;  // Should be 10000 SM3/day in SI
    std::cout << "FIELD ORAT target: " << metric_rate(field_orat_target) << " SM3/day" << std::endl;

    // Get reduction rates
    const auto& field_red = gs.production_reduction_rates("FIELD");
    const auto& plat_red = gs.production_reduction_rates("PLAT");
    std::cout << "FIELD reduction[oil]: " << metric_rate(field_red[oil_phase_pos]) << " SM3/day" << std::endl;
    std::cout << "PLAT reduction[oil]: " << metric_rate(plat_red[oil_phase_pos]) << " SM3/day" << std::endl;

    // Compute expected guide rate fractions
    // MANI and MANI2 are siblings under PLAT
    // Guide rates: WELL-A=5000, WELL-B=4000
    const double mani_guide = 5000.0;
    const double mani2_guide = 4000.0;
    const double total_guide = mani_guide + mani2_guide;
    const double mani_fraction = mani_guide / total_guide;
    std::cout << "Guide rate fraction for MANI: " << mani_fraction << std::endl;

    // Expected calculation (simplified, may not match actual algorithm)
    const double accumulated_eff = E_MANI * E_PLAT;
    std::cout << "Accumulated efficiency: " << accumulated_eff << std::endl;

    // If target = (FIELD_target - FIELD_reduction + addback) * fraction / eff_factor
    // This is approximate - actual algorithm is more complex

    std::cout << "\nNote: The actual algorithm involves multiple iterations through" << std::endl;
    std::cout << "the group chain with local_reduction_level logic. Use GDB to" << std::endl;
    std::cout << "trace the exact calculation in checkGroupConstraintsProd()." << std::endl;

    std::cout << "\n=== Isolated Test Complete ===" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
