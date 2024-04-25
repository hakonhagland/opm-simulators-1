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
#ifndef OPM_PY_BLACKOIL_SIMULATOR_IMPL_HEADER_INCLUDED
#define OPM_PY_BLACKOIL_SIMULATOR_IMPL_HEADER_INCLUDED
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMain.hpp>
//#include <opm/simulators/flow/python/PyFluidState.hpp>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#include <opm/simulators/flow/python/PyBlackOilSimulator.hpp>
// NOTE: EXIT_SUCCESS, EXIT_FAILURE is defined in cstdlib
#include <cstdlib>
#include <stdexcept>
#include <string>


namespace py = pybind11;

namespace Opm::Pybind {
template <class Version>
PyBlackOilSimulator<Version>::PyBlackOilSimulator( const std::string &deck_filename)
    : deck_filename_{deck_filename}
{
}
template <class Version>
PyBlackOilSimulator<Version>::PyBlackOilSimulator(
    std::shared_ptr<Opm::Deck> deck,
    std::shared_ptr<Opm::EclipseState> state,
    std::shared_ptr<Opm::Schedule> schedule,
    std::shared_ptr<Opm::SummaryConfig> summary_config
)
    : deck_{std::move(deck)}
    , eclipse_state_{std::move(state)}
    , schedule_{std::move(schedule)}
    , summary_config_{std::move(summary_config)}
{
}

// Public methods alphabetically sorted
// ------------------------------------

template <class Version>
void PyBlackOilSimulator<Version>::advance(int report_step)
{
    while (currentStep() < report_step) {
        step();
    }
}

template <class Version>
bool PyBlackOilSimulator<Version>::checkSimulationFinished()
{
    return getFlowMain().getSimTimer()->done();
}

// This returns the report step number that will be executed next time step()
//   is called.
template <class Version>
int PyBlackOilSimulator<Version>::currentStep()
{
    return getFlowMain().getSimTimer()->currentStepNum();
    // NOTE: this->simulator_->episodeIndex() would also return the current
    // report step number, but this number is always delayed by 1 step relative
    // to this->flow_main_->getSimTimer()->currentStepNum()
    // See details in runStep() in file SimulatorFullyImplicitBlackoilEbos.hpp
}

template <class Version>
py::array_t<double> PyBlackOilSimulator<Version>::getCellVolumes() {
    auto vector = getMaterialState().getCellVolumes();
    return py::array(vector.size(), vector.data());
}

template <class Version>
double PyBlackOilSimulator<Version>::getDT() {
    return getFlowMain().getPreviousReportStepSize();
}

template <class Version>
py::array_t<double> PyBlackOilSimulator<Version>::getPorosity()
{
    auto vector = getMaterialState().getPorosity();
    return py::array(vector.size(), vector.data());
}

template <class Version>
py::array_t<double>
PyBlackOilSimulator<Version>::
getFluidStateVariable(const std::string &name) const
{
    auto vector = getFluidState().getFluidStateVariable(name);
    return py::array(vector.size(), vector.data());
}

template <class Version>
py::array_t<double>
PyBlackOilSimulator<Version>::
getPrimaryVariable(const std::string &variable) const
{
    auto vector = getFluidState().getPrimaryVariable(variable);
    return py::array(vector.size(), vector.data());
}

template <class Version>
py::array_t<int>
PyBlackOilSimulator<Version>::
getPrimaryVarMeaning(const std::string &variable) const
{
    auto vector = getFluidState().getPrimaryVarMeaning(variable);
    return py::array(vector.size(), vector.data());
}

template <class Version>
std::map<std::string, int>
PyBlackOilSimulator<Version>::
getPrimaryVarMeaningMap(const std::string &variable) const
{

    return getFluidState().getPrimaryVarMeaningMap(variable);
}

template <class Version>
int PyBlackOilSimulator<Version>::run()
{
    auto main_object = Opm::Main( this->deck_filename_ );
    return main_object.runStatic<Opm::Properties::TTag::FlowProblemTPFA>();
}

template <class Version>
void PyBlackOilSimulator<Version>::setPorosity( py::array_t<double,
    py::array::c_style | py::array::forcecast> array)
{
    std::size_t size_ = array.size();
    const double *poro = array.data();
    getMaterialState().setPorosity(poro, size_);
}

template <class Version>
void
PyBlackOilSimulator<Version>::
setPrimaryVariable(
    const std::string &variable,
    py::array_t<double,
    py::array::c_style | py::array::forcecast> array
)
{
    std::size_t size_ = array.size();
    const double *data = array.data();
    getFluidState().setPrimaryVariable(variable, data, size_);
}

template <class Version>
int PyBlackOilSimulator<Version>::step()
{
    if (!this->has_run_init_) {
        throw std::logic_error("step() called before step_init()");
    }
    if (this->has_run_cleanup_) {
        throw std::logic_error("step() called after step_cleanup()");
    }
    if(checkSimulationFinished()) {
        throw std::logic_error("step() called, but simulation is done");
    }
    //if (this->debug_)
    //    this->mainEbos_->getSimTimer()->report(std::cout);
    auto result = getFlowMain().executeStep();
    return result;
}

template <class Version>
int PyBlackOilSimulator<Version>::stepCleanup()
{
    this->has_run_cleanup_ = true;
    return getFlowMain().executeStepsCleanup();
}

template <class Version>
int PyBlackOilSimulator<Version>::stepInit()
{

    if (this->has_run_init_) {
        // Running step_init() multiple times is not implemented yet,
        if (this->has_run_cleanup_) {
            throw std::logic_error("step_init() called again");
        }
        else {
            return EXIT_SUCCESS;
        }
    }
    if (this->deck_) {
        this->main_ = std::make_unique<Opm::Main>(
            this->deck_->getDataFile(),
            this->eclipse_state_,
            this->schedule_,
            this->summary_config_
        );
    }
    else {
        this->main_ = std::make_unique<Opm::Main>( this->deck_filename_ );
    }
    int exit_code = EXIT_SUCCESS;
    this->flow_main_ = this->main_->initFlowBlackoil(exit_code);
    if (this->flow_main_) {
        int result = this->flow_main_->executeInitStep();
        this->has_run_init_ = true;
        this->simulator_ = this->flow_main_->getSimulatorPtr();
        this->fluid_state_ = std::make_unique<PyFluidState<TypeTag>>(this->simulator_);
        this->material_state_ = std::make_unique<PyMaterialState<TypeTag>>(this->simulator_);
        return result;
    }
    else {
        return exit_code;
    }
}

// Private methods alphabetically sorted
// ------------------------------------

template <class Version>
Opm::FlowMain<typename PyBlackOilSimulator<Version>::TypeTag>&
         PyBlackOilSimulator<Version>::getFlowMain() const
{
    if (this->flow_main_) {
        return *this->flow_main_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

template <class Version>
PyFluidState<typename PyBlackOilSimulator<Version>::TypeTag>&
PyBlackOilSimulator<Version>::
getFluidState() const
{
    if (this->fluid_state_) {
        return *this->fluid_state_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMainEbos object" );
    }
}

template <class Version>
PyMaterialState<typename PyBlackOilSimulator<Version>::TypeTag>&
PyBlackOilSimulator<Version>::getMaterialState() const
{
    if (this->material_state_) {
        return *this->material_state_;
    }
    else {
        throw std::runtime_error("BlackOilSimulator not initialized: "
            "Cannot get reference to FlowMain object" );
    }
}

} // namespace Opm::Pybind


#endif // OPM_PY_BLACKOIL_SIMULATOR_IMPL_HEADER_INCLUDED
