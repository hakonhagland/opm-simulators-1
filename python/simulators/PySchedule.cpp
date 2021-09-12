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

#include "config.h"
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

namespace Opm::Pybind
{
const Opm::Well& PySchedule::getWell(
    const std::string& name, const size_t& timestep ) {
    if (auto wrapper = this->schedule_wrapper_.lock()) {
        auto &schedule = wrapper->getSchedule();
        try {
            return schedule.getWell( name, timestep );
        } catch( const std::invalid_argument& e ) {
            throw py::key_error( name );
        }
    }
    else {
        throw std::runtime_error("Cannot get well: BlackOilSimulator does not exist");
    }
}

void export_PySchedule(py::module& m)
{
    using namespace Opm::Pybind;
    py::class_<PySchedule>(m, "PySchedule")
        .def(py::init<std::weak_ptr<PyScheduleWrapper>>())
        .def( "get_well", &PySchedule::getWell);
}
} //namespace Opm::Pybind

