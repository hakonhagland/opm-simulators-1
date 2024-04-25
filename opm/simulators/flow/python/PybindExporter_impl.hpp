/*
  Copyright 2024 Equinor ASA.

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

#ifndef OPM_PYBIND11_EXPORTER_IMPL_HEADER_INCLUDED
#define OPM_PYBIND11_EXPORTER_IMPL_HEADER_INCLUDED

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMain.hpp>
//#include <opm/simulators/flow/python/PyFluidState.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>

#include <pybind11/pybind11.h>
#include <opm/simulators/flow/python/PybindExporter.hpp>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>
// NOTE: This file will be generated at compile time and placed in the build directory
// See python/generate_docstring_hpp.py, and python/simulators/CMakeLists.txt for details
#include <PyBlackOilSimulatorDoc.hpp>

namespace Opm::Pybind {

template <class Version>
void export_all(py::module& m, const std::string& class_name) {
    export_PyBlackOilSimulator<Version>(m, class_name);
}

// Exported functions
template <class Version>
void export_PyBlackOilSimulator(py::module& m, const std::string& class_name)
{
    using namespace Opm::Pybind::DocStrings;

    py::class_<PyBlackOilSimulator<Version>>(m, class_name.c_str())
        .def(py::init<const std::string&>(),
             PyBlackOilSimulator_filename_constructor_docstring)
        .def(py::init<
             std::shared_ptr<Opm::Deck>,
             std::shared_ptr<Opm::EclipseState>,
             std::shared_ptr<Opm::Schedule>,
             std::shared_ptr<Opm::SummaryConfig>>(),
             PyBlackOilSimulator_objects_constructor_docstring)
        .def("advance", &PyBlackOilSimulator<Version>::advance, advance_docstring, py::arg("report_step"))
        .def("check_simulation_finished", &PyBlackOilSimulator<Version>::checkSimulationFinished,
             checkSimulationFinished_docstring)
        .def("current_step", &PyBlackOilSimulator<Version>::currentStep, currentStep_docstring)
        .def("get_cell_volumes", &PyBlackOilSimulator<Version>::getCellVolumes, getCellVolumes_docstring)
        .def("get_dt", &PyBlackOilSimulator<Version>::getDT, getDT_docstring)
        .def("get_fluidstate_variable", &PyBlackOilSimulator<Version>::getFluidStateVariable,
            py::return_value_policy::copy, py::arg("name"))
        .def("get_porosity", &PyBlackOilSimulator<Version>::getPorosity, getPorosity_docstring)
        .def("get_primary_variable_meaning", &PyBlackOilSimulator<Version>::getPrimaryVarMeaning,
            py::return_value_policy::copy, py::arg("variable"))
        .def("get_primary_variable_meaning_map", &PyBlackOilSimulator<Version>::getPrimaryVarMeaningMap,
            py::return_value_policy::copy, py::arg("variable"))
        .def("get_primary_variable", &PyBlackOilSimulator<Version>::getPrimaryVariable,
            py::return_value_policy::copy, py::arg("variable"))
        .def("run", &PyBlackOilSimulator<Version>::run, run_docstring)
        .def("set_porosity", &PyBlackOilSimulator<Version>::setPorosity, setPorosity_docstring, py::arg("array"))
        .def("set_primary_variable", &PyBlackOilSimulator<Version>::setPrimaryVariable,
            py::arg("variable"), py::arg("value"))
        .def("step", &PyBlackOilSimulator<Version>::step, step_docstring)
        .def("step_cleanup", &PyBlackOilSimulator<Version>::stepCleanup, stepCleanup_docstring)
        .def("step_init", &PyBlackOilSimulator<Version>::stepInit, stepInit_docstring);
}

} // namespace Opm::Pybind
#endif //OPM_PYBIND11_EXPORTER_IMPL_HEADER_INCLUDED
