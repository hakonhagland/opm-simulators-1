/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <config.h>

#include <ebos/ebos.hh>
#include <ebos/eclgenericvanguard.hh>

#include <opm/common/utility/Serializer.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>

#include <opm/models/blackoil/blackoilprimaryvariables.hh>

#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/utils/SerializationPackers.hpp>

#define BOOST_TEST_MODULE TestRestartSerialization
#define BOOST_TEST_NO_MAIN

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/test/unit_test.hpp>

template<class T>
std::tuple<T,int,int> PackUnpack(T& in)
{
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(in);
    const size_t pos1 = ser.position();
    T out{};
    ser.unpack(out);
    const size_t pos2 = ser.position();

    return std::make_tuple(std::move(out), pos1, pos2);
}

#define TEST_FOR_TYPE_NAMED_OBJ(TYPE, NAME, OBJ) \
BOOST_AUTO_TEST_CASE(NAME) \
{ \
    auto val1 = Opm::TYPE::OBJ(); \
    auto val2 = PackUnpack(val1); \
    BOOST_CHECK_MESSAGE(std::get<1>(val2) == std::get<2>(val2), "Packed size differ from unpack size for " #TYPE); \
    BOOST_CHECK_MESSAGE(val1 == std::get<0>(val2), "Deserialized " #TYPE " differ"); \
}

#define TEST_FOR_TYPE_NAMED(TYPE, NAME) \
    TEST_FOR_TYPE_NAMED_OBJ(TYPE, NAME, serializationTestObject)

#define TEST_FOR_TYPE(TYPE) \
    TEST_FOR_TYPE_NAMED(TYPE, TYPE)

TEST_FOR_TYPE(HardcodedTimeStepControl)
TEST_FOR_TYPE(PIDAndIterationCountTimeStepControl)
TEST_FOR_TYPE(PIDTimeStepControl)
TEST_FOR_TYPE(SimpleIterationCountTimeStepControl)
TEST_FOR_TYPE(SimulatorTimer)

namespace Opm { using ATE = AdaptiveTimeSteppingEbos<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosHardcoded, serializationTestObjectHardcoded)
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosPID, serializationTestObjectPID)
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosPIDIt, serializationTestObjectPIDIt)
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosSimple, serializationTestObjectSimple)

namespace Opm { using BPV = BlackOilPrimaryVariables<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE_NAMED(BPV, BlackoilPrimaryVariables)

BOOST_AUTO_TEST_CASE(EclGenericVanguard)
{
    auto in_params = Opm::EclGenericVanguard::serializationTestParams();
    Opm::EclGenericVanguard val1(std::move(in_params));
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(val1);
    const size_t pos1 = ser.position();
    Opm::EclGenericVanguard::SimulationModelParams out_params;
    out_params.setupTime_ = 0.0;
    out_params.actionState_ = std::make_unique<Opm::Action::State>();
    out_params.udqState_ = std::make_unique<Opm::UDQState>();
    out_params.eclSchedule_ = std::make_shared<Opm::Schedule>();
    out_params.summaryState_ = std::make_unique<Opm::SummaryState>();
    Opm::EclGenericVanguard val2(std::move(out_params));
    ser.unpack(val2);
    const size_t pos2 = ser.position();

    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for EclGenericVanguard");
    BOOST_CHECK_MESSAGE(val1 == val2, "Deserialized EclGenericVanguard differ");
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}