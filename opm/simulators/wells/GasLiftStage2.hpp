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

#ifndef OPM_GASLIFT_STAGE2_HEADER_INCLUDED
#define OPM_GASLIFT_STAGE2_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
// NOTE: BlackoilWellModel.hpp includes ourself (GasLiftStage2.hpp), so we need
//   to forward declare BlackoilWellModel for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class BlackoilWellModel;
}
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftStage2 {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using BlackoilWellModel = Opm::BlackoilWellModel<TypeTag>;
        using MPICommunicator = Dune::MPIHelper::MPICommunicator;
        using CollectiveCommunication = Dune::CollectiveCommunication<MPICommunicator>;
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
    public:
        GasLiftStage2(
            const BlackoilWellModel &well_model,
            const Simulator &ebos_simulator,
            DeferredLogger &deferred_logger,
            const WellState &well_state
        );
        void runOptimize();
    private:
        void optimizeGroup(const Opm::Group &group);

        DeferredLogger &deferred_logger_;
        const Simulator &ebos_simulator_;
        const BlackoilWellModel &well_model_;
        const WellState &well_state_;

        int report_step_idx_;
        const SummaryState &summary_state_;
        const CollectiveCommunication &comm_;
        const Schedule &schedule_;
        const PhaseUsage &phase_usage_;
    };

} // namespace Opm

#include "GasLiftStage2_impl.hpp"

#endif // OPM_GASLIFT_STAGE2_HEADER_INCLUDED
