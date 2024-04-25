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

#ifndef OPM_PYBIND11_EXPORTER_HEADER_INCLUDED
#define OPM_PYBIND11_EXPORTER_HEADER_INCLUDED

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
//#include <pybind11/embed.h>

namespace py = pybind11;

namespace Opm::Pybind {
    template <class Version>
    void export_all(py::module& m, const std::string& class_name);
    template <class Version>
    void export_PyBlackOilSimulator(py::module& m, const std::string& class_name);
}

#ifndef OPM_PYBIND11_EXPORTER_IMPL_HEADER_INCLUDED
#include "PybindExporter_impl.hpp"
#endif

#endif //OPM_PYBIND11_EXPORTER_HEADER_INCLUDED
