#include "config.h"

#include <pybind11/pybind11.h>
#include <opm/simulators/flow/python/PybindExporter.hpp>
#include <opm/simulators/flow/FlowMain.hpp>

namespace Opm {

std::unique_ptr<FlowMain<Properties::TTag::FlowProblemTPFA>>
flowBlackoilTpfaMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    return std::make_unique<FlowMain<Properties::TTag::FlowProblemTPFA>>(
        argc, argv, outputCout, outputFiles);
}

}

struct Version1 {};
PYBIND11_MODULE(simulators, m)
{
    Opm::Pybind::export_all<Version1>(m, "BlackOilSimulator");
}
