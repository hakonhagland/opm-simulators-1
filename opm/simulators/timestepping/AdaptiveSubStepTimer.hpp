/*
  Copyright 2024 Equinor AS

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

#ifndef OPM_ADAPTIVE_SUB_STEP_TIMER_HEADER_INCLUDED
#define OPM_ADAPTIVE_SUB_STEP_TIMER_HEADER_INCLUDED

#include <cassert>
#include <iosfwd>
#include <vector>
#include <limits>
#include <algorithm>
#include <memory>
#include <numeric>

namespace Opm
{

    /////////////////////////////////////////////////////////////
    ///
    /// \brief Simulation timer for adaptive substep timestepping
    ///
    //////////////////////////////////////////////////////////////
    class AdaptiveSubStepTimer
    {
    public:
        /// \brief constructor
        ///  \param simulation_time_elapsed  time elapsed since start of simulation
        ///  \param current_step_length      length of current time step
        ///  \param last_step_taken            last suggested time step
        ///  \param max_time_step              maximum time step allowed
        AdaptiveSubStepTimer( double simulation_time_elapsed,
                              double current_step_length,
                              const double last_step_taken,
                              const double max_time_step = std::numeric_limits<double>::max() );

        /// \brief advance time by currentStepLength
        AdaptiveSubStepTimer& operator++ ();

        /// \brief advance time by currentStepLength
        void advance() { this->operator++ (); }

        /// \brief provide and estimate for new time step size
        void provideTimeStepEstimate( const double dt_estimate );

        /// \brief \copydoc SimulationTimer::currentStepNum
        int currentStepNum () const;

        /// \brief \copydoc SimulationTimer::currentStepLength
        double currentStepLength () const;

        // \brief Set next step length
        void setCurrentStepLength(double dt);

        /// \brief \copydoc SimulationTimer::totalTime
        double totalTime() const;

        /// \brief \copydoc SimulationTimer::simulationTimeElapsed
        double simulationTimeElapsed() const;

        /// \brief \copydoc SimulationTimer::done
        bool done () const;

        /// \brief return average step length used so far
        double averageStepLength() const;

        /// \brief return max step length used so far
        double maxStepLength () const;

        /// \brief return min step length used so far
        double minStepLength () const;

        /// \brief Previous step length. This is the length of the step that
        ///        was taken to arrive at this time.
        double stepLengthTaken () const;

        /// \brief report start and end time as well as used steps so far
        void report(std::ostream& os) const;

        /// \brief Return true if last time step failed
        bool lastStepFailed() const {return last_step_failed_;}

        /// \brief tell the timestepper whether timestep failed or not
        void setLastStepFailed(bool last_step_failed) {last_step_failed_ = last_step_failed;}

    protected:
        const double start_time_;
        const double total_time_;
        const double max_time_step_;

        double current_time_;
        double dt_;
        int current_step_;

        std::vector< double > steps_;
        bool last_step_failed_;

    };

} // namespace Opm

#endif // OPM_SUB_STEP_TIMER_HEADER_INCLUDED
