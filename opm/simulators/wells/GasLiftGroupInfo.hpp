/*
  Copyright 2021 Equinor ASA.

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

#ifndef OPM_GASLIFT_GROUP_INFO_HEADER_INCLUDED
#define OPM_GASLIFT_GROUP_INFO_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
// NOTE: BlackoilWellModel.hpp includes ourself (GasLiftStage2.hpp), so we need
//   to forward declare BlackoilWellModel for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class BlackoilWellModel;
}
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftGroupInfo
    {
        class GroupRates;
        using BlackoilWellModel = Opm::BlackoilWellModel<TypeTag>;
        using WellInterfacePtr = typename BlackoilWellModel::WellInterfacePtr;
        using WellInterfaceRawPtr = typename BlackoilWellModel::WellInterfaceRawPtr;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        // NOTE: In the Well2GroupMap below, in the std::vector value we store
        //    pairs of group names and efficiency factors. The efficiency factors
        //    are the product of the wells efficiency factor and all the efficiency
        //    factors of the child groups of the group all the way down
        //    to the well group.
        using Well2GroupMap =
            std::map<std::string, std::vector<std::pair<std::string,double>>>;
        using GroupRateMap =
            std::map<std::string, GroupRates>;
        // TODO: same definition with WellInterface, and
        //  WellStateFullyImplicitBlackoil eventually they should go
        //  to a common header file.
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
    public:
        GasLiftGroupInfo(
            const BlackoilWellModel &well_model,
            const Simulator &ebos_simulator,
            DeferredLogger &deferred_logger,
            WellState &well_state);
        std::vector<std::pair<std::string,double>> &getWellGroups(
            const std::string &well_name);
        void initialize();
        double alqRate(const std::string &group_name);
        double gasRate(const std::string &group_name);
        std::optional<double> gasTarget(const std::string &group_name);
        std::optional<double> maxAlq(const std::string &group_name);
        double oilRate(const std::string &group_name);
        std::optional<double> oilTarget(const std::string &group_name);
        void update(const std::string &well_name,
            double delta_oil, double delta_gas, double delta_alq);
        const Well2GroupMap &wellGroupMap() { return well_group_map_; }
    private:
        bool checkDoGasLiftOptimization_(const std::string &well_name);
        bool checkNewtonIterationIdxOk_(const std::string &well_name);
        void displayDebugMessage_(const std::string &msg);
        void displayDebugMessage_(const std::string &msg, const std::string &well_name);
        std::pair<double, double> getProducerWellRates_(WellInterfaceRawPtr well);
        std::tuple<double, double, double> initializeGroupRatesRecursive_(
            const Group &group);
        void initializeWell2GroupMapRecursive_(
            const Group &group, std::vector<std::string> &group_names,
            std::vector<double> &group_efficiency, double cur_efficiency);
        class GroupRates {
        public:
            GroupRates( double oil_rate, double gas_rate, double alq,
                std::optional<double> oil_target,
                std::optional<double> gas_target,
                std::optional<double> total_gas,
                std::optional<double> max_alq
            ) :
                oil_rate_{oil_rate},
                gas_rate_{gas_rate},
                alq_{alq},
                oil_target_{oil_target},
                gas_target_{gas_target},
                total_gas_{total_gas},
                max_alq_{max_alq}
            {}
            double alq() const { return alq_; }
            double gasRate() const { return gas_rate_; }
            std::optional<double> gasTarget() const { return gas_target_; }
            std::optional<double> maxAlq() const { return max_alq_; }
            double oilRate() const { return oil_rate_; }
            std::optional<double> oilTarget() const { return oil_target_; }
            void update(double delta_oil, double delta_gas, double delta_alq)
            {
                oil_rate_ += delta_oil;
                gas_rate_ += delta_gas;
                alq_ += delta_alq;
            }
        private:
            double oil_rate_;
            double gas_rate_;
            double alq_;
            std::optional<double> oil_target_;
            std::optional<double> gas_target_;
            std::optional<double> total_gas_;
            std::optional<double> max_alq_;
        };

        const BlackoilWellModel &well_model_;
        const Simulator &ebos_simulator_;
        DeferredLogger &deferred_logger_;
        WellState &well_state_;
        const Schedule &schedule_;
        const SummaryState &summary_state_;
        int report_step_idx_;
        const GasLiftOpt& glo_;
        GroupRateMap group_rate_map_;
        Well2GroupMap well_group_map_;
        bool debug;

        // Optimize only wells under THP control
        bool optimize_only_thp_wells_ = true;

    };

#include "GasLiftGroupInfo_impl.hpp"

} // namespace Opm

#endif // OPM_GASLIFT_GROUP_INFO_INCLUDED
