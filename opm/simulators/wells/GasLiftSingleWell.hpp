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

#ifndef OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
#define OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
// NOTE: StandardWell.hpp includes ourself (GasLiftSingleWell.hpp), so we need
//   to forward declare StandardWell for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class StandardWell;
}
#include <opm/simulators/wells/StandardWell.hpp>

#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftSingleWell
    {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using StdWell = Opm::StandardWell<TypeTag>;
        // TODO: same definition with WellInterface, and
        //  WellStateFullyImplicitBlackoil eventually they should go
        //  to a common header file.
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
        struct OptimizeState;
        class Stage2State;
    public:
        GasLiftSingleWell(
            const StdWell &std_well,
            const Simulator &ebos_simulator,
            const SummaryState &summary_state,
            DeferredLogger &deferred_logger,
            WellState &well_state
        );
        struct GradInfo;
        void addOrSubtractAlqIncrement(const GradInfo &gi, bool increase);
        std::optional<GradInfo> calcIncOrDecGradient(bool increase) const;
        void runOptimize();
        const std::string& name() const {return well_name_; }

        struct GradInfo
        {
            GradInfo() {}
            GradInfo(double grad_, double new_oil_rate_, bool oil_is_limited_,
                     double new_gas_rate_, bool gas_is_limited_,
                     double alq_, bool alq_is_limited_) :
                grad{grad_},
                new_oil_rate{new_oil_rate_},
                oil_is_limited{oil_is_limited_},
                new_gas_rate{new_gas_rate_},
                gas_is_limited{gas_is_limited_},
                alq{alq_},
                alq_is_limited{alq_is_limited_} {}
            double grad;
            double new_oil_rate;
            bool oil_is_limited;
            double new_gas_rate;
            bool gas_is_limited;
            double alq;
            bool alq_is_limited;
        };

    private:
        std::pair<double, bool> addOrSubtractAlqIncrement_(
            double alq, bool increase) const;
        double calcEcoGradient_(double oil_rate, double new_oil_rate,
            double gas_rate, double new_gas_rate) const;
        bool checkWellRatesViolated_(
            std::vector<double> &potentials,
            const std::function<bool(double, double, const std::string &)> &callback);
        std::optional<double> computeBhpAtThpLimit_(double alq) const;
        void computeInitialWellRates_();
        void computeWellRates_(double bhp, std::vector<double> &potentials) const;
        void debugCheckNegativeGradient_(double grad, double alq, double new_alq,
            double oil_rate, double new_oil_rate, double gas_rate,
            double new_gas_rate, bool increase) const;
        void debugShowIterationInfo_(OptimizeState &state, double alq);
        void debugShowStartIteration_(double alq, bool increase);
        void displayDebugMessage_(const std::string &msg) const;
        void displayWarning_(std::string warning);
        std::pair<double, bool> getBhpWithLimit_(double bhp) const;
        std::pair<double, bool> getGasRateWithLimit_(
            const std::vector<double> &potentials) const;
        std::pair<double, bool> getOilRateWithLimit_(
            const std::vector<double> &potentials) const;
        void logSuccess_(double alq);
        std::optional<double> runOptimizeLoop_(bool increase);
        void setAlqMaxRate_(const GasLiftOpt::Well &well);
        void setAlqMinRate_(const GasLiftOpt::Well &well);
        bool stage2CheckRateAlreadyLimited_(bool increase) const;
        std::optional<double> tryIncreaseLiftGas_();
        std::optional<double> tryDecreaseLiftGas_();
        void updateWellStateAlqFixedValue_(const GasLiftOpt::Well &well);
        bool useFixedAlq_(const GasLiftOpt::Well &well);
        void warnMaxIterationsExceeded_();

        DeferredLogger &deferred_logger_;
        const Simulator &ebos_simulator_;
        std::vector<double> potentials_;
        const StdWell &std_well_;
        const SummaryState &summary_state_;
        WellState &well_state_;
        std::string well_name_;
        const Well &ecl_well_;
        const Well::ProductionControls controls_;
        bool debug;  // extra debug output

        double alpha_w_;
        double alpha_g_;
        double eco_grad_;
        int gas_pos_;
        bool has_run_init_ = false;
        double increment_;
        double max_alq_;
        int max_iterations_;
        double min_alq_;
        int oil_pos_;
        bool optimize_;
        double orig_alq_;
        int water_pos_;
        std::optional<Stage2State> stage2_state_;

        struct OptimizeState
        {
            OptimizeState( GasLiftSingleWell &parent_, bool increase_ ) :
                parent{parent_},
                increase{increase_},
                it{0},
                stop_iteration{false},
                bhp{-1}
            {}

            GasLiftSingleWell &parent;
            bool increase;
            int it;
            bool stop_iteration;
            double bhp;

            std::pair<double,bool> addOrSubtractAlqIncrement(double alq);
            bool checkAlqOutsideLimits(double alq, double oil_rate);
            bool checkEcoGradient(double gradient);
            bool checkOilRateExceedsTarget(double oil_rate);
            bool checkRate(double rate, double limit, const std::string &rate_str) const;
            bool checkWellRatesViolated(std::vector<double> &potentials);
            bool computeBhpAtThpLimit(double alq);
            double getBhpWithLimit();
            void warn_(std::string msg) {parent.displayWarning_(msg);}
        };

        class Stage2State
        {
        public:
            Stage2State(double oil_rate, bool oil_is_limited,
                        double gas_rate, bool gas_is_limited,
                        double alq, bool alq_is_limited,
                        bool increase) :
                oil_rate_{oil_rate},
                oil_is_limited_{oil_is_limited},
                gas_rate_{gas_rate},
                gas_is_limited_{gas_is_limited},
                alq_{alq},
                alq_is_limited_{alq_is_limited},
                increase_{increase} {}
            double alq() const { return alq_; }
            bool alqIsLimited() const { return alq_is_limited_; }
            bool gasIsLimited() const { return gas_is_limited_; }
            double gasRate() const { return gas_rate_; }
            bool oilIsLimited() const { return oil_is_limited_; }
            double oilRate() const { return oil_rate_; }
            bool increase() const { return increase_; }
            void update(double oil_rate, bool oil_is_limited,
                        double gas_rate, bool gas_is_limited,
                        double alq, bool alq_is_limited,
                        bool increase)
            {
                oil_rate_ = oil_rate;
                oil_is_limited_ = oil_is_limited;
                gas_rate_ = gas_rate;
                gas_is_limited_ = gas_is_limited;
                alq_ = alq;
                alq_is_limited_ = alq_is_limited;
                increase_ = increase;
            }
        private:
            double oil_rate_;
            bool oil_is_limited_;
            double gas_rate_;
            bool gas_is_limited_;
            double alq_;
            bool alq_is_limited_;
            bool increase_;
        };
    };

#include "GasLiftSingleWell_impl.hpp"

} // namespace Opm


#endif // OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
