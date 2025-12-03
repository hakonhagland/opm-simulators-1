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


#ifndef OPM_TARGETCALCULATOR_HEADER_INCLUDED
#define OPM_TARGETCALCULATOR_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>

#include <optional>
#include <string>
#include <vector>

namespace Opm {

class DeferredLogger;
template<class Scalar> class GroupState;
template<typename IndexTraits> class PhaseUsageInfo;
template<typename Scalar, typename IndexTraits> class GroupStateHelper;

namespace GroupStateHelpers
 {

/// Based on a group control mode, extract or calculate rates, and
/// provide other conveniences.
template<typename Scalar, typename IndexTraits>
class TargetCalculator
{
public:
    using GroupStateHelperType = Opm::GroupStateHelper<Scalar, IndexTraits>;

    TargetCalculator(const GroupStateHelperType& groupStateHelper,
                     const std::vector<Scalar>& resv_coeff,
                     const Group& group);

    template <typename RateType>
    RateType calcModeRateFromRates(const std::vector<RateType>& rates,
                                   const bool rates_are_reservoir = false) const
    {
      return calcModeRateFromRates(rates.data(), rates_are_reservoir);
    }

    template <typename RateType>
    RateType calcModeRateFromRates(const RateType* rates,
                                   const bool rates_are_reservoir = false) const;

    Scalar groupTarget(DeferredLogger& deferred_logger) const;

    GuideRateModel::Target guideTargetMode() const;

    /// @brief Check if the control mode is RESV (reservoir volume rate)
    /// @return true if control mode is RESV, false otherwise
    bool isResvMode() const { return this->cmode_ == Group::ProductionCMode::RESV; }

private:
    Group::ProductionCMode cmode_;
    const GroupStateHelperType& groupStateHelper_;
    const std::vector<Scalar>& resv_coeff_;
    const Group& group_;
};

/// Based on a group control mode, extract or calculate rates, and
/// provide other conveniences.
template<typename Scalar, typename IndexTraits>
class InjectionTargetCalculator
{
public:
    using GroupStateHelperType = GroupStateHelper<Scalar, IndexTraits>;

    InjectionTargetCalculator(const GroupStateHelperType& groupStateHelper,
                              const std::vector<Scalar>& resv_coeff,
                              const Group& group,
                              const Phase& injection_phase,
                              DeferredLogger& deferred_logger);

    template <typename RateVec>
    auto calcModeRateFromRates(const RateVec& rates) const
    {
        return rates[pos_];
    }

    Scalar groupTarget(DeferredLogger& deferred_logger) const;

    GuideRateModel::Target guideTargetMode() const;

private:
    const GroupStateHelperType& groupStateHelper_;
    const std::vector<Scalar>& resv_coeff_;
    const Group& group_;
    const Phase& injection_phase_;
    Group::InjectionCMode cmode_;
    int pos_;
    GuideRateModel::Target target_;
};

} // namespace GroupStateHelpers


} // namespace Opm

#endif
