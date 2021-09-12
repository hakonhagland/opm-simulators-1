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

#ifndef OPM_PY_BLACKOIL_SIMULATOR_HEADER_INCLUDED
#define OPM_PY_BLACKOIL_SIMULATOR_HEADER_INCLUDED

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/flow/python/Pybind11Exporter.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <opm/simulators/flow/python/PySchedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

namespace Opm::Pybind {
class PyBlackOilSimulator
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblem;
private:
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

public:
    PyBlackOilSimulator( const std::string &deckFilename);
    py::array_t<double> getPorosity();
    PySchedule *getSchedule();
    int run();
    void setPorosity(
         py::array_t<double, py::array::c_style | py::array::forcecast> array);
    int step();
    int stepInit();
    int stepCleanup();
    const Opm::FlowMainEbos<TypeTag>& getFlowMainEbos() const;

private:
    const std::string deckFilename_;
    bool hasRunInit_ = false;
    bool hasRunCleanup_ = false;

    // This *must* be declared before other pointers
    // to simulator objects. This in order to deinitialize
    // MPI at the correct time (ie after the other objects).
    std::unique_ptr<Opm::Main> main_;

    std::unique_ptr<Opm::FlowMainEbos<TypeTag>> mainEbos_;
    Simulator *ebosSimulator_;
    std::unique_ptr<PyMaterialState<TypeTag>> materialState_;
    std::shared_ptr<PyScheduleWrapper> scheduleWrapper_;
};

class PyScheduleWrapper {
public:
    PyScheduleWrapper(const PyBlackOilSimulator& sim) : sim_(sim)
    {
    }
    const Opm::Schedule& getSchedule();
private:
    const PyBlackOilSimulator& sim_;
};

} // namespace Opm::Pybind
#endif // OPM_PY_BLACKOIL_SIMULATOR_HEADER_INCLUDED
