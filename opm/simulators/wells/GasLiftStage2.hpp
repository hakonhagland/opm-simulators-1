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
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
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
        using GasLiftSingleWell = Opm::GasLiftSingleWell<TypeTag>;
        using GLiftOptWells = typename BlackoilWellModel::GLiftOptWells;
        using GLiftProdWells = typename BlackoilWellModel::GLiftProdWells;
        using GradPair = std::pair<std::string, double>;
        using GradInfo = typename GasLiftSingleWell::GradInfo;
        using GradMap = std::map<std::string, GradInfo>;
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
    public:
        GasLiftStage2(
            const BlackoilWellModel &well_model,
            const Simulator &ebos_simulator,
            DeferredLogger &deferred_logger,
            const WellState &well_state,
            GLiftProdWells &prod_wells,
            GLiftOptWells &glift_wells
        );
        void runOptimize();
    private:
        std::optional<double> calcIncOrDecGrad_(
            const std::string name, GasLiftSingleWell &gs_well, bool increase);
        void displayDebugMessage_(const std::string &msg);
        void displayDebugMessage_(const std::string &msg, const std::string &group_name);
        void displayWarning_(const std::string &msg, const std::string &group_name);
        std::vector<GasLiftSingleWell *> getGroupGliftWells_(
            const Opm::Group &group);
        void getGroupGliftWellsRecursive_(
            const Opm::Group &group, std::vector<GasLiftSingleWell *> &wells);
        void optimizeGroup_(const Opm::Group &group);
        void optimizeGroupsRecursive_(const Opm::Group &group);
        void saveGrad_(GradMap &map, const std::string name, GradInfo &grad);
        void saveDecGrad_(const std::string name, GradInfo &grad);
        void saveIncGrad_(const std::string name, GradInfo &grad);

        DeferredLogger &deferred_logger_;
        const Simulator &ebos_simulator_;
        const BlackoilWellModel &well_model_;
        const WellState &well_state_;
        GLiftProdWells &prod_wells_;
        GLiftOptWells &stage1_wells_;

        int report_step_idx_;
        const SummaryState &summary_state_;
        const CollectiveCommunication &comm_;
        const Schedule &schedule_;
        const PhaseUsage &phase_usage_;
        const GasLiftOpt& glo_;
        GradMap inc_grads_;
        GradMap dec_grads_;
        bool debug_;
        int max_iterations_ = 1000;

        struct OptimizeState {
            OptimizeState( GasLiftStage2 &parent_, const Opm::Group &group_ ) :
                parent{parent_},
                group{group_},
                it{0}
            {}
            GasLiftStage2 &parent;
            const Opm::Group &group;
            int it;

            using GradInfo = typename GasLiftStage2::GradInfo;
            using GradPair = typename GasLiftStage2::GradPair;
            using GradMap = typename GasLiftStage2::GradMap;
            std::pair<std::optional<GradPair>,std::optional<GradPair>>
              calculateEcoGradients(
                 std::vector<GasLiftSingleWell *> &wells);
            bool checkAtLeastTwoWells(std::vector<GasLiftSingleWell *> &wells);
            void debugShowIterationInfo();
            void redistributeALQ( const std::string well1, double grad1,
                const std::string well2, double grad2);

        private:
            void addOrRemoveALQincrement_(
                GradMap &grad_map, const std::string well_name, bool add);
            std::optional<GradInfo> calcIncOrDecGrad_(
                const std::string well_name,
                const GasLiftSingleWell &gs_well, bool increase);
            void displayDebugMessage_(const std::string &msg);
            void displayWarning_(const std::string &msg);
            void sortGradients_(std::vector<GradPair> &grads);

        };
    };

#include "GasLiftStage2_impl.hpp"

} // namespace Opm

#endif // OPM_GASLIFT_STAGE2_HEADER_INCLUDED
