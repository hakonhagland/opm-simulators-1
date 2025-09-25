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

#include <opm/simulators/flow/rescoup/SendSlaveGroupData.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

namespace Opm {

template <class Scalar, class IndexTraits>
SendSlaveGroupData<Scalar, IndexTraits>::
SendSlaveGroupData(GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
                   const GroupState<Scalar>& group_state,
                   int report_step)
    : guide_rate_handler_(guide_rate_handler)
    , group_state_(group_state)
    , report_step_(report_step)
{
}

template <class Scalar, class IndexTraits>
void
SendSlaveGroupData<Scalar, IndexTraits>::
sendGroupRatesAndPotentialsToMaster()
{
    SlaveGroupRates slave_rates;
    slave_rates.report_step = report_step_;

    // TODO: Determine group name from reservoir coupling configuration
    // For now, use placeholder
    slave_rates.group_name = "SLAVE_GROUP";

    // Gather all required data
    gatherPotentials_(slave_rates);
    gatherProductionRates_(slave_rates);
    gatherInjectionRates_(slave_rates);

    // Send to master
    sendToMaster_(slave_rates);
}

template <class Scalar, class IndexTraits>
void
SendSlaveGroupData<Scalar, IndexTraits>::
gatherPotentials_(SlaveGroupRates& slave_rates)
{
    // TODO: Extract potentials gathering logic from GuideRateHandler
    // This replaces the logic currently in GuideRateHandler::sendSlaveGroupPotentialsToMaster()

    // Placeholder implementation - needs to be moved from GuideRateHandler
    slave_rates.potentials.rate[static_cast<std::size_t>(ReservoirCoupling::Potentials<Scalar>::Phase::Oil)] = 0.0;
    slave_rates.potentials.rate[static_cast<std::size_t>(ReservoirCoupling::Potentials<Scalar>::Phase::Gas)] = 0.0;
    slave_rates.potentials.rate[static_cast<std::size_t>(ReservoirCoupling::Potentials<Scalar>::Phase::Water)] = 0.0;
}

template <class Scalar, class IndexTraits>
void
SendSlaveGroupData<Scalar, IndexTraits>::
gatherProductionRates_(SlaveGroupRates& slave_rates)
{
    // Issue 5 & 6: Production rates needed for guide rate fractions and REIN calculations
    // TODO: Extract production rates from GroupState for slave groups

    const auto& phase_usage = guide_rate_handler_.phaseUsage();
    const int num_phases = phase_usage.numActivePhases();

    slave_rates.production_rates.resize(num_phases, 0.0);

    // TODO: Get actual group name from reservoir coupling configuration
    // and extract production rates from group_state_.production_rates(group_name)

    // Placeholder - needs actual group name and rate extraction
    for (int phase = 0; phase < num_phases; ++phase) {
        slave_rates.production_rates[phase] = 0.0;
    }
}

template <class Scalar, class IndexTraits>
void
SendSlaveGroupData<Scalar, IndexTraits>::
gatherInjectionRates_(SlaveGroupRates& slave_rates)
{
    // Issue 4 & 7: Injection rates needed for guide rates and GPMAINT
    const auto& phase_usage = guide_rate_handler_.phaseUsage();
    const int num_phases = phase_usage.numActivePhases();

    slave_rates.injection_surface_rates.resize(num_phases, 0.0);
    slave_rates.injection_reservoir_rates.resize(num_phases, 0.0);
    slave_rates.injection_vrep_rate = 0.0;

    // TODO: Get actual group name from reservoir coupling configuration
    // and extract rates from:
    // - group_state_.injection_surface_rates(group_name)
    // - group_state_.injection_reservoir_rates(group_name)
    // - group_state_.injection_vrep_rate(group_name)

    // Placeholder - needs actual group name and rate extraction
    for (int phase = 0; phase < num_phases; ++phase) {
        slave_rates.injection_surface_rates[phase] = 0.0;
        slave_rates.injection_reservoir_rates[phase] = 0.0;
    }
}

template <class Scalar, class IndexTraits>
void
SendSlaveGroupData<Scalar, IndexTraits>::
sendToMaster_(const SlaveGroupRates& slave_rates)
{
    // TODO: Implement MPI communication to master
    // This should:
    // 1. Get master communicator from reservoir coupling slave
    // 2. Serialize SlaveGroupRates using MPI traits
    // 3. Send with appropriate message tag

    OpmLog::debug("SendSlaveGroupData: Sending slave group data to master");
    OpmLog::debug("  Group: " + slave_rates.group_name);
    OpmLog::debug("  Report step: " + std::to_string(slave_rates.report_step));
    OpmLog::debug("  Production phases: " + std::to_string(slave_rates.production_rates.size()));
    OpmLog::debug("  Injection phases: " + std::to_string(slave_rates.injection_surface_rates.size()));

    // Placeholder - actual MPI sending logic needed
}

// Explicit template instantiations
template class SendSlaveGroupData<double, BlackOilDefaultFluidSystemIndices>;
template class SendSlaveGroupData<float, BlackOilDefaultFluidSystemIndices>;

} // namespace Opm