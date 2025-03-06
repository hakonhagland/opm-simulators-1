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

#ifndef OPM_WELL_MODEL_RESERVOIR_COUPLING_HANDLER_HPP
#define OPM_WELL_MODEL_RESERVOIR_COUPLING_HANDLER_HPP

#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm {
template <typename TypeTag> class BlackoilWellModel;

template<typename TypeTag>
class WellModelReservoirCouplingHandler {
public:
    WellModelReservoirCouplingHandler(
        BlackoilWellModel<TypeTag>& well_model, DeferredLogger& deferred_logger);

    void sendSlaveGroupPotentialsToMaster();
private:
    BlackoilWellModel<TypeTag>& well_model_;
    DeferredLogger& deferred_logger_;
};

} // namespace Opm

#include "WellModelReservoirCouplingHandler_impl.hpp"

#endif // OPM_WELL_MODEL_RESERVOIR_COUPLING_HANDLER_HPP
