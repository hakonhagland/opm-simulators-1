/*
  Copyright 2014 IRIS AS

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
#ifndef OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED
#define OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED

#include <cassert>
#include <iosfwd>
#include <vector>
#include <limits>
#include <algorithm>
#include <memory>
#include <numeric>

#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/simulators/timestepping/AdaptiveSubStepTimer.hpp>

namespace Opm
{

    /////////////////////////////////////////////////////////
    ///
    /// \brief Simulation timer for adaptive time stepping
    ///
    /////////////////////////////////////////////////////////
    class AdaptiveSimulatorTimer : public SimulatorTimerInterface, public AdaptiveSubStepTimer
    {
    public:
        /// \brief constructor
        ///  \param start_date               start date of simulation
        ///  \param simulation_time_elapsed  time elapsed since start of simulation
        ///  \param current_step_length      length of current time step
        ///  \param report_step_num          current report step
        ///  \param lastStepTaken            last suggested time step
        ///  \param maxTimeStep              maximum time step allowed
        AdaptiveSimulatorTimer( boost::posix_time::ptime start_date,
                                double simulation_time_elapsed,
                                double current_step_length,
                                int report_step_num,
                                const double last_step_taken,
                                const double max_time_step = std::numeric_limits<double>::max() );

        // ----------------------------------------------------------------------
        // Implementation of pure virtual functions from SimulatorTimerInterface
        // ----------------------------------------------------------------------

        /// \brief advance time by currentStepLength
        void advance() override { AdaptiveSubStepTimer::advance(); }

        /// return copy of object
        virtual std::unique_ptr< SimulatorTimerInterface > clone() const override;

        /// \brief \copydoc SimulationTimer::currentStepLength
        double currentStepLength () const override { return AdaptiveSubStepTimer::currentStepLength(); }

        /// \brief \copydoc SimulationTimer::currentStepNum
        int currentStepNum () const override { return AdaptiveSubStepTimer::currentStepNum(); }

        /// \brief \copydoc SimulationTimer::done
        bool done () const override { return AdaptiveSubStepTimer::done(); }

        /// \brief Whether this is the first step
        bool initialStep () const override;

        /// \brief Return true if last time step failed
        bool lastStepFailed() const override {return AdaptiveSubStepTimer::lastStepFailed();}

        /// \brief \copydoc SimulationTimer::simulationTimeElapsed
        double simulationTimeElapsed() const override { return AdaptiveSubStepTimer::simulationTimeElapsed(); }

        /// \brief Previous step length. This is the length of the step that
        ///        was taken to arrive at this time.
        double stepLengthTaken () const override { return AdaptiveSubStepTimer::stepLengthTaken(); }

        /// \brief start date time of simulation
        boost::posix_time::ptime startDateTime() const override;

        // --------------------------------------------------------------------------------------
        // Overriding the default implementation of virtual functions in SimulatorTimerInterface
        // --------------------------------------------------------------------------------------

        /// \brief return current report step
        int reportStepNum() const override { return report_step_; }

        // -----------------------------------------------------
        // Methods not part of the SimulatorTimerInterface
        // -----------------------------------------------------


        double averageStepLength() const { return AdaptiveSubStepTimer::averageStepLength(); }

        /// \brief return average step length used so far
        /// \brief return max step length used so far
        double maxStepLength () const { return AdaptiveSubStepTimer::maxStepLength(); }

        /// \brief return min step length used so far
        double minStepLength () const { return AdaptiveSubStepTimer::minStepLength(); }

        /// \brief provide and estimate for new time step size
        void provideTimeStepEstimate( const double dt_estimate ) {
            AdaptiveSubStepTimer::provideTimeStepEstimate( dt_estimate ); 
        }

        /// \brief report start and end time as well as used steps so far
        void report(std::ostream& os) const {
            AdaptiveSubStepTimer::report(os);
        }

        // \brief Set next step length
        void setCurrentStepLength(double dt) {
            AdaptiveSubStepTimer::setCurrentStepLength(dt);
        }

        /// \brief \copydoc SimulationTimer::totalTime
        double totalTime() const { return AdaptiveSubStepTimer::totalTime(); }

    protected:
        std::shared_ptr<boost::posix_time::ptime> start_date_time_;
        const int report_step_;
    };

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
