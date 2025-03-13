/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.
  Copyright 2023 Inria

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
#ifndef OPM_PYBASEMAIN_HEADER_INCLUDED
#define OPM_PYBASEMAIN_HEADER_INCLUDED

#include <opm/simulators/flow/Main.hpp>

#include <string>

namespace Opm {
class PyBaseMain : public Main
{
public:
    void setArguments(const std::vector<std::string>& args)
    {
        if (args.empty()) {
            return;
        }

        // We have the two arguments previously setup (binary name and input
        // case name) by the main class plus whichever args are in the
        // parameter that was passed from the python side.
        this->argc_ = 2 + args.size();

        // Setup our vector of char*'s
        argv_python_.resize(2 + args.size());
        argv_python_[0] = argv_[0];
        argv_python_[1] = argv_[1];
        for (std::size_t i = 0; i < args.size(); ++i) {
            argv_python_[i+2] = const_cast<char*>(args[i].c_str());
        }

        // Finally set the main class' argv pointer to the combined
        // parameter list.
        this->argv_ = argv_python_.data();
    }

protected:
    std::vector<char*> argv_python_{};
};  // class PyBaseMain
}  // namespace Opm

#endif