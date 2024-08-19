/*
  Copyright (c) 2014 IRIS AS

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ostream>
#include <numeric>
#include <vector>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/AdaptiveSubStepTimer.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace Opm
{
    AdaptiveSimulatorTimer::
    AdaptiveSimulatorTimer( boost::posix_time::ptime start_date,
                            double simulation_time_elapsed,
                            double current_step_length,
                            int report_step_num,
                            const double last_step_taken,
                            const double max_time_step )
        : start_date_time_{std::make_shared<boost::posix_time::ptime>(start_date)}
        , AdaptiveSubStepTimer{simulation_time_elapsed, current_step_length, last_step_taken, max_time_step}
        , report_step_{report_step_num}
    {
    }

    bool AdaptiveSimulatorTimer::initialStep () const
    {
        return ( reportStepNum() == 0 ) && ( currentStepNum() == 0 );
    }

    boost::posix_time::ptime AdaptiveSimulatorTimer::startDateTime() const
    {
        return *start_date_time_;
    }

    /// return copy of object
    std::unique_ptr< SimulatorTimerInterface >
    AdaptiveSimulatorTimer::clone() const
    {
        return std::make_unique<AdaptiveSimulatorTimer>(*this);
    }

} // namespace Opm
