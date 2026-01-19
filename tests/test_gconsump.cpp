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
/*!\
 * \file
 * \brief Comprehensive test suite for GCONSUMP efficiency factor handling
 *
 * This test suite verifies hierarchical efficiency factor application in GCONSUMP calculations.
 *
 * \section test_hierarchy Group Hierarchy
 * \code
 *                       FIELD (GEFAC=1.0)
 *                         |
 *                      PLAT-A (GEFAC=0.9)
 *                         |
 *                      SUB-B (GEFAC=0.8, GCONSUMP consumption=100, import=50)
 *                         |
 *                      WELL-C
 * \endcode
 *
 * \section test_objectives Primary Test Objectives
 * - Verify GCONSUMP rates are properly adjusted by hierarchical efficiency factors
 * - Confirm efficiency factors are applied multiplicatively (Eclipse specification)
 * - Validate update_gconsump() correctly accumulates efficiency-adjusted rates
 * - Test unit conversion between SI and metric rates
 * - Demonstrate proper tolerance levels for different assertion types
 *
 * \section expected_behavior Expected Behavior (Following REIN Pattern)
 * - Raw rates: consumption=100, import=50 (from Eclipse deck)
 * - SUB-B rates: consumption=100, import=50 (raw, not affected by own efficiency)
 * - PLAT-A rates: consumption=80, import=40 (SUB-B rates × SUB-B efficiency 0.8)
 * - FIELD rates: consumption=72, import=36 (PLAT-A rates × PLAT-A efficiency 0.9)
 *
 * \section test_fixture Test Fixture Features
 * - ToleranceAndUnitFixture provides:
 *   - Multiple tolerance levels (tight, rate, algorithm)
 *   - Unit conversion helpers (SI ↔ metric)
 *   - Assertion helpers with appropriate tolerances
 */

#define BOOST_TEST_MODULE GConsumpEfficiencyFactor
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/simulators/wells/GroupState.hpp>

#include <boost/test/unit_test.hpp>

using namespace Opm;

//! \brief Test fixture with unit conversion helpers and tolerance checks
struct ToleranceAndUnitFixture
{
    //! Relative tolerance for exact mathematical relationships
    static constexpr double tight_tol = 1e-12;

    //! Relative tolerance for rate values with minor numerical differences
    static constexpr double rate_tol = 1e-8;

    //! Relative tolerance for algorithm results (e.g., scale factors)
    static constexpr double algo_tol = 0.1;

    //! Convert SI rate (m³/s) to metric rate (SM3/day)
    static double metric_rate(double si_rate)
    {
        using namespace Opm::unit;
        return convert::to(si_rate, cubic(meter) / day);
    }

    //! Convert metric rate (SM3/day) to SI rate (m³/s)
    static double si_rate(double metric_rate)
    {
        using namespace Opm::unit;
        return convert::from(metric_rate, cubic(meter) / day);
    }

    //! Check that two values are close within tight tolerance (for exact relationships)
    static void checkClose(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, tight_tol);
    }

    //! Check that two rate values are close
    static void checkRate(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, rate_tol);
    }

    //! Check that two algorithm results are close (looser tolerance)
    static void checkAlgo(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, algo_tol);
    }
};

BOOST_FIXTURE_TEST_SUITE(GConsumpEfficiencyFactorTests, ToleranceAndUnitFixture)

//! Test that demonstrates the GCONSUMP efficiency factor correct implementation
//!
//! This test verifies that update_gconsump() applies efficiency factors following the REIN pattern.
//! Expected behavior (following REIN efficiency pattern):
//! - SUB-B has GEFAC=0.8, PLAT-A has GEFAC=0.9
//! - SUB-B GCONSUMP: consumption=100, import=50
//! - Expected SUB-B rates: consumption=100, import=50 (raw rates, not affected by own efficiency)
//! - Expected PLAT-A rates: consumption=80, import=40 (SUB-B rates × SUB-B efficiency)
//! - Expected FIELD rates: consumption=72, import=36 (PLAT-A rates × PLAT-A efficiency)
//!
//! This follows the same pattern as REIN where groups are not affected by their own efficiency,
//! but parents apply child efficiency factors when accumulating child contributions.
BOOST_AUTO_TEST_CASE(TestGconsumpEfficiencyFactorBug)
{
    // Load the test deck with GCONSUMP and GEFAC specification
    const std::string deck_file = "GCONSUMP.DATA";
    Parser parser;
    auto deck = parser.parseFile(deck_file);
    EclipseState eclipseState(deck);
    Schedule schedule(deck, eclipseState);
    const int report_step = 0;
    const auto& sched_state = schedule[report_step];
    const auto& gconsump = sched_state.gconsump();

    const std::size_t num_phases = 3; // oil, water, gas

    // Initialize group state
    GroupState<double> group_state(num_phases);

    // Create summary state for UDQ support
    SummaryState summary_state(schedule.getStartTime());

    // Update GCONSUMP rates - this is where the bug manifests
    group_state.update_gconsump(schedule, report_step, summary_state);

    const auto& group_sub_b = schedule.getGroup("SUB-B", report_step);
    const auto& group_plat_a = schedule.getGroup("PLAT-A", report_step);
    const auto& group_field = schedule.getGroup("FIELD", report_step);

    const auto& sub_b_efficiency = group_sub_b.getGroupEfficiencyFactor();
    const auto& plat_a_efficiency = group_plat_a.getGroupEfficiencyFactor();
    const auto& field_efficiency = group_field.getGroupEfficiencyFactor();

    // Validate efficiency factors are as expected
    checkClose(sub_b_efficiency, 0.8);
    checkClose(plat_a_efficiency, 0.9);
    checkClose(field_efficiency, 1.0);

    // Note: Cumulative efficiency calculation is for reference only
    // The actual efficiency application follows the REIN pattern:
    // - Groups are NOT affected by their own efficiency
    // - Parent groups apply child efficiency when accumulating child contributions
    const double cumulative_efficiency = sub_b_efficiency * plat_a_efficiency;
    checkClose(cumulative_efficiency, 0.72);

    // Get the efficiency-adjusted rates from update_gconsump()
    const auto& sub_b_rates = group_state.gconsump_rates("SUB-B");
    const auto& plat_a_rates = group_state.gconsump_rates("PLAT-A");
    const auto& field_rates = group_state.gconsump_rates("FIELD");

    // Raw GCONSUMP values for SUB-B from deck
    const auto& group_gc = gconsump.get("SUB-B", summary_state);
    const double raw_consumption = metric_rate(group_gc.consumption_rate);
    const double raw_import = metric_rate(group_gc.import_rate);

    // Eclipse-compliant expected values (correct efficiency pattern)
    // SUB-B: raw rates (not affected by its own efficiency)
    const double expected_sub_b_consumption = raw_consumption;    // 100.0
    const double expected_sub_b_import = raw_import;              // 50.0

    // PLAT-A: SUB-B's rates × SUB-B's efficiency (0.8)
    const double expected_plat_a_consumption = raw_consumption * sub_b_efficiency;  // 80.0
    const double expected_plat_a_import = raw_import * sub_b_efficiency;            // 40.0

    // FIELD: PLAT-A's rates × PLAT-A's efficiency (0.9)
    const double expected_field_consumption = expected_plat_a_consumption * plat_a_efficiency;  // 72.0
    const double expected_field_import = expected_plat_a_import * plat_a_efficiency;            // 36.0

    // Convert actual SI rates to metric for direct comparison with expected metric values
    const double actual_sub_b_consumption_metric = metric_rate(sub_b_rates.first);
    const double actual_sub_b_import_metric = metric_rate(sub_b_rates.second);
    const double actual_plat_a_consumption_metric = metric_rate(plat_a_rates.first);
    const double actual_plat_a_import_metric = metric_rate(plat_a_rates.second);
    const double actual_field_consumption_metric = metric_rate(field_rates.first);
    const double actual_field_import_metric = metric_rate(field_rates.second);

    // Check if the ratios are correct (efficiency factors might affect units but should preserve ratios)
    const double actual_ratio = (actual_sub_b_import_metric > 0) ? actual_sub_b_consumption_metric / actual_sub_b_import_metric : 0;
    const double expected_ratio = raw_consumption / raw_import; // 2.0

    // Validate that the consumption/import ratio is preserved (efficiency should not change relative proportions)
    checkClose(actual_ratio, expected_ratio);

    std::cout << "\n=== GCONSUMP Efficiency Factor Test Results ===" << std::endl;
    std::cout << "Expected (Eclipse-compliant with correct efficiency pattern) [SM3/day]:" << std::endl;
    std::cout << "  SUB-B: consumption=" << expected_sub_b_consumption << ", import=" << expected_sub_b_import << std::endl;
    std::cout << "  PLAT-A: consumption=" << expected_plat_a_consumption << ", import=" << expected_plat_a_import << std::endl;
    std::cout << "  FIELD: consumption=" << expected_field_consumption << ", import=" << expected_field_import << std::endl;

    std::cout << "\nActual (current implementation) [SM3/day]:" << std::endl;
    std::cout << "  SUB-B: consumption=" << actual_sub_b_consumption_metric << ", import=" << actual_sub_b_import_metric << std::endl;
    std::cout << "  PLAT-A: consumption=" << actual_plat_a_consumption_metric << ", import=" << actual_plat_a_import_metric << std::endl;
    std::cout << "  FIELD: consumption=" << actual_field_consumption_metric << ", import=" << actual_field_import_metric << std::endl;

    std::cout << "\n=== Efficiency Factor Validation ===" << std::endl;
    std::cout << "SUB-B efficiency factor: " << sub_b_efficiency << std::endl;
    std::cout << "PLAT-A efficiency factor: " << plat_a_efficiency << std::endl;
    std::cout << "FIELD efficiency factor: " << field_efficiency << std::endl;
    std::cout << "Cumulative efficiency (SUB-B): " << cumulative_efficiency << std::endl;

    std::cout << "\n=== Diagnostic Information ===" << std::endl;
    std::cout << "Raw GCONSUMP specification [SM3/day]: consumption=" << raw_consumption << ", import=" << raw_import << std::endl;
    std::cout << "Expected efficiency pattern [SM3/day]:" << std::endl;
    std::cout << "  SUB-B (raw): " << expected_sub_b_consumption << " / " << expected_sub_b_import << std::endl;
    std::cout << "  PLAT-A (×0.8): " << expected_plat_a_consumption << " / " << expected_plat_a_import << std::endl;
    std::cout << "  FIELD (×0.9): " << expected_field_consumption << " / " << expected_field_import << std::endl;
    std::cout << "Expected ratio (consumption/import): " << expected_ratio << std::endl;
    std::cout << "Actual ratio (consumption/import): " << actual_ratio << std::endl;
    std::cout << "Efficiency factors correctly applied: " << (std::abs(actual_ratio - expected_ratio) < tight_tol ? "YES" : "NO") << std::endl;

    // Note: All comparisons are done in metric units (SM3/day) for direct comparison with Eclipse deck values
    // The efficiency factors follow the REIN pattern where groups are not affected by their own efficiency

    std::cout << "\n=== Testing SUB-B rates (where GCONSUMP is specified) ===" << std::endl;

    // SUB-B rates should be raw (NOT affected by its own efficiency)
    checkRate(actual_sub_b_consumption_metric, expected_sub_b_consumption);
    checkRate(actual_sub_b_import_metric, expected_sub_b_import);

    std::cout << "\n=== Testing PLAT-A rates (hierarchical accumulation) ===" << std::endl;

    // PLAT-A accumulates SUB-B's contribution with SUB-B's efficiency factor
    checkRate(actual_plat_a_consumption_metric, expected_plat_a_consumption);
    checkRate(actual_plat_a_import_metric, expected_plat_a_import);

    std::cout << "\n=== Testing FIELD rates (top-level accumulation) ===" << std::endl;

    // FIELD accumulates PLAT-A's contribution with PLAT-A's efficiency factor
    checkRate(actual_field_consumption_metric, expected_field_consumption);
    checkRate(actual_field_import_metric, expected_field_import);
}

//! Test for groups without GCONSUMP - should return zero rates
BOOST_AUTO_TEST_CASE(TestGconsumpZeroRates)
{
    // Load the test deck
    const std::string deck_file = "GCONSUMP.DATA";
    Parser parser;
    auto deck = parser.parseFile(deck_file);
    EclipseState eclipseState(deck);
    Schedule schedule(deck, eclipseState);

    const int report_step = 0;
    const std::size_t num_phases = 3; // oil, water, gas

    // Initialize group state
    GroupState<double> group_state(num_phases);
    SummaryState summary_state(schedule.getStartTime());

    // Update GCONSUMP rates
    group_state.update_gconsump(schedule, report_step, summary_state);

    // WELL-C has no GCONSUMP specified - should return zero rates
    const auto& well_c_rates = group_state.gconsump_rates("WELL-C");

    BOOST_CHECK_EQUAL(well_c_rates.first, 0.0);
    BOOST_CHECK_EQUAL(well_c_rates.second, 0.0);
}

//! Test complex GCONSUMP efficiency factor hierarchy
//!
//! This test validates the implementation with a realistic, complex group structure:
//! - Multiple sibling groups
//! - GCONSUMP at multiple hierarchy levels
//! - Mixed scenarios (groups with/without GCONSUMP)
//! - Complex efficiency factor accumulation
//!
//! Structure:
//! FIELD (GEFAC=1.0, GCONSUMP=20/10)
//! ├── PLAT-A (GEFAC=0.9, GCONSUMP=30/15)
//! │   ├── SUB-B (GEFAC=0.8, GCONSUMP=100/50)
//! │   │   └── WELL-C
//! │   └── SUB-D (GEFAC=0.7, no GCONSUMP)
//! │       └── WELL-E
//! └── PLAT-F (GEFAC=0.85, no GCONSUMP)
//!     └── SUB-G (GEFAC=0.75, GCONSUMP=80/40)
//!         └── WELL-H
BOOST_AUTO_TEST_CASE(TestGconsumpComplexHierarchy)
{
    // Load the complex test deck
    const std::string deck_file = "GCONSUMP_COMPLEX.DATA";
    Parser parser;
    auto deck = parser.parseFile(deck_file);
    EclipseState eclipseState(deck);
    Schedule schedule(deck, eclipseState);

    const int report_step = 0;
    const std::size_t num_phases = 3; // oil, water, gas

    // Initialize group state
    GroupState<double> group_state(num_phases);
    SummaryState summary_state(schedule.getStartTime());

    // Update GCONSUMP rates
    group_state.update_gconsump(schedule, report_step, summary_state);

    // Get all group rates
    const auto& sub_b_rates = group_state.gconsump_rates("SUB-B");
    const auto& sub_d_rates = group_state.gconsump_rates("SUB-D");
    const auto& plat_a_rates = group_state.gconsump_rates("PLAT-A");
    const auto& field_rates = group_state.gconsump_rates("FIELD");

    // Convert to metric for comparison
    const auto sub_b_consumption_metric = metric_rate(sub_b_rates.first);
    const auto sub_b_import_metric = metric_rate(sub_b_rates.second);
    const auto sub_d_consumption_metric = metric_rate(sub_d_rates.first);
    const auto sub_d_import_metric = metric_rate(sub_d_rates.second);
    const auto plat_a_consumption_metric = metric_rate(plat_a_rates.first);
    const auto plat_a_import_metric = metric_rate(plat_a_rates.second);
    const auto field_consumption_metric = metric_rate(field_rates.first);
    const auto field_import_metric = metric_rate(field_rates.second);

    // Expected values (following REIN efficiency pattern)
    const double expected_sub_b_consumption = 100.0;    // Own GCONSUMP only
    const double expected_sub_b_import = 50.0;

    const double expected_sub_d_consumption = 0.0;      // No GCONSUMP
    const double expected_sub_d_import = 0.0;

    // PLAT-A: own 30/15 + children (SUB-B: 100×0.8, SUB-D: 0×0.7) = 30+80+0/15+40+0 = 110/55
    const double expected_plat_a_consumption = 30.0 + (100.0 * 0.8) + (0.0 * 0.7);  // 110.0
    const double expected_plat_a_import = 15.0 + (50.0 * 0.8) + (0.0 * 0.7);        // 55.0

    // FIELD: own 20/10 + children (PLAT-A: 110×0.9) = 20+99/10+49.5 = 119/59.5
    const double expected_field_consumption = 20.0 + (110.0 * 0.9);  // 119.0
    const double expected_field_import = 10.0 + (55.0 * 0.9);        // 59.5

    std::cout << "\n=== Complex GCONSUMP Hierarchy Test Results ===" << std::endl;
    std::cout << "Expected (Eclipse-compliant) [SM3/day]:" << std::endl;
    std::cout << "  SUB-B: consumption=" << expected_sub_b_consumption << ", import=" << expected_sub_b_import << std::endl;
    std::cout << "  SUB-D: consumption=" << expected_sub_d_consumption << ", import=" << expected_sub_d_import << std::endl;
    std::cout << "  PLAT-A: consumption=" << expected_plat_a_consumption << ", import=" << expected_plat_a_import << std::endl;
    std::cout << "  FIELD: consumption=" << expected_field_consumption << ", import=" << expected_field_import << std::endl;

    std::cout << "\nActual (implementation) [SM3/day]:" << std::endl;
    std::cout << "  SUB-B: consumption=" << sub_b_consumption_metric << ", import=" << sub_b_import_metric << std::endl;
    std::cout << "  SUB-D: consumption=" << sub_d_consumption_metric << ", import=" << sub_d_import_metric << std::endl;
    std::cout << "  PLAT-A: consumption=" << plat_a_consumption_metric << ", import=" << plat_a_import_metric << std::endl;
    std::cout << "  FIELD: consumption=" << field_consumption_metric << ", import=" << field_import_metric << std::endl;

    // Test leaf group with own GCONSUMP
    std::cout << "\n=== Testing leaf group with GCONSUMP ===" << std::endl;
    checkRate(sub_b_consumption_metric, expected_sub_b_consumption);
    checkRate(sub_b_import_metric, expected_sub_b_import);

    // Test group without GCONSUMP
    std::cout << "\n=== Testing group without GCONSUMP ===" << std::endl;
    checkRate(sub_d_consumption_metric, expected_sub_d_consumption);
    checkRate(sub_d_import_metric, expected_sub_d_import);

    // Test intermediate group (with both own GCONSUMP and children)
    std::cout << "\n=== Testing intermediate group with own+children GCONSUMP ===" << std::endl;
    checkRate(plat_a_consumption_metric, expected_plat_a_consumption);
    checkRate(plat_a_import_metric, expected_plat_a_import);

    // Test field-level accumulation (most complex)
    std::cout << "\n=== Testing field-level accumulation ===" << std::endl;
    checkRate(field_consumption_metric, expected_field_consumption);
    checkRate(field_import_metric, expected_field_import);
}

BOOST_AUTO_TEST_SUITE_END()
