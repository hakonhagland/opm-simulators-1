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

#ifndef OPM_PY_WELL_HEADER_INCLUDED
#define OPM_PY_WELL_HEADER_INCLUDED

#include <opm/models/utils/propertysystem.hh>

#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

namespace Opm::Pybind
{
    class PySimulatorWrapper;
    class PyWell {
    public:
        PyWell(std::weak_ptr<PySimulatorWrapper> simulator_wrapper,
            const Opm::Well& well)
            :
               simulator_wrapper_(simulator_wrapper),
               well_{well}
        {
        }
        const std::string& getName();

    private:
        std::weak_ptr<PySimulatorWrapper> simulator_wrapper_;
        const Opm::Well &well_;
    };

}
#endif // OPM_PY_WELL_HEADER_INCLUDED
