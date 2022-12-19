/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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


#ifndef OPM_MULTISEGMENTWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_PRIMARY_VARIABLES_HEADER_INCLUDED

#include <opm/material/densead/Evaluation.hpp>

#include <array>
#include <vector>

namespace Opm
{

template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class MultisegmentWellPrimaryVariables
{
public:
    // TODO: for now, not considering the polymer, solvent and so on to simplify the development process.

    // TODO: we need to have order for the primary variables and also the order for the well equations.
    // sometimes, they are similar, while sometimes, they can have very different forms.

    // Table showing the primary variable indices, depending on what phases are present:
    //
    //         WOG     OG     WG     WO    W/O/G (single phase)
    // WQTotal   0      0      0      0                       0
    // WFrac     1  -1000      1      1                   -1000
    // GFrac     2      1  -1000  -1000                   -1000
    // Spres     3      2      2      2                       1

    static constexpr bool has_water = (Indices::waterSwitchIdx >= 0);
    static constexpr bool has_gas = (Indices::compositionSwitchIdx >= 0);
    static constexpr bool has_oil = (Indices::numPhases - has_gas - has_water) > 0;

    // In the implementation, one should use has_wfrac_variable
    // rather than has_water to check if you should do something
    // with the variable at the WFrac location, similar for GFrac.
    static constexpr bool has_wfrac_variable = has_water && Indices::numPhases > 1;
    static constexpr bool has_gfrac_variable = has_gas && has_oil;

    static constexpr int WQTotal = 0;
    static constexpr int WFrac = has_wfrac_variable ? 1 : -1000;
    static constexpr int GFrac = has_gfrac_variable ? has_wfrac_variable + 1 : -1000;
    static constexpr int SPres = has_wfrac_variable + has_gfrac_variable + 1;

    //  the number of well equations  TODO: it should have a more general strategy for it
    static constexpr int numWellEq = Indices::numPhases + 1;

    using EvalWell = DenseAd::Evaluation<double, /*size=*/Indices::numEq + numWellEq>;

    MultisegmentWellPrimaryVariables(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well)
        : well_(well)
    {}

    //! \brief Resize values and evaluations.
    void resize(const int numSegments);

    //! \brief Initialize evaluations from values.
    void init();

    //! \brief Copy values from well state.
    void update(const WellState& well_state);

    // the values for the primary varibles
    // based on different solutioin strategies, the wells can have different primary variables
    std::vector<std::array<double, numWellEq> > value_;

    // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
    std::vector<std::array<EvalWell, numWellEq> > evaluation_;

    //! \brief Handle non-reasonable fractions due to numerical overshoot.
    void processFractions(const int seg);

private:
    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_; //!< Reference to well interface
};

}

#endif // OPM_MULTISEGMENTWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
