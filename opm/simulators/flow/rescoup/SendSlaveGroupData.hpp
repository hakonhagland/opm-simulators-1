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

#ifndef OPM_SEND_SLAVE_GROUP_DATA_HPP
#define OPM_SEND_SLAVE_GROUP_DATA_HPP

#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/wells/GroupState.hpp>

namespace Opm {

template <class Scalar, class IndexTraits> class GuideRateHandler;

/// \brief Handles sending comprehensive slave group data to master process
///
/// This class consolidates the sending of potentials and rates from slave
/// groups to the master process. It replaces the scattered logic in
/// GuideRateHandler and provides a clean, single-responsibility interface
/// for slave data transmission.
///
/// The class gathers:
/// - Group potentials (existing functionality)
/// - Production rates (for guide rate fractions and REIN calculations)
/// - Injection surface rates (for GPMAINT pressure maintenance)
/// - Injection reservoir rates (for guide rate calculations)
/// - Injection vrep rates (for guide rate NETV targets)
template <class Scalar, class IndexTraits>
class SendSlaveGroupData {
public:
    using SlaveGroupRates = ReservoirCoupling::SlaveGroupRates<Scalar>;

    /// \brief Constructor
    /// \param guide_rate_handler Reference to guide rate handler
    /// \param group_state Current group state with rate data
    /// \param report_step Current report step index
    SendSlaveGroupData(GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
                       const GroupState<Scalar>& group_state,
                       int report_step);

    /// \brief Send comprehensive group data to master process
    ///
    /// This method:
    /// 1. Gathers potentials from guide rate handler
    /// 2. Collects production rates from group state
    /// 3. Collects injection rates (surface, reservoir, vrep) from group state
    /// 4. Packages all data into SlaveGroupRates structure
    /// 5. Sends via MPI to master process
    void sendGroupRatesAndPotentialsToMaster();

private:
    /// \brief Gather potentials for all slave groups
    /// \param slave_rates Output structure to populate with potentials
    void gatherPotentials_(SlaveGroupRates& slave_rates);

    /// \brief Gather production rates for all slave groups
    /// \param slave_rates Output structure to populate with production rates
    void gatherProductionRates_(SlaveGroupRates& slave_rates);

    /// \brief Gather injection rates for all slave groups
    /// \param slave_rates Output structure to populate with injection rates
    void gatherInjectionRates_(SlaveGroupRates& slave_rates);

    /// \brief Send packed data to master process via MPI
    /// \param slave_rates Complete data structure to send
    void sendToMaster_(const SlaveGroupRates& slave_rates);

    // Member variables
    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler_;
    const GroupState<Scalar>& group_state_;
    int report_step_;
};

} // namespace Opm

#endif // OPM_SEND_SLAVE_GROUP_DATA_HPP