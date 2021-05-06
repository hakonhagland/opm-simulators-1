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

template<typename TypeTag>
GasLiftGroupInfo<TypeTag>::
GasLiftGroupInfo(
    const BlackoilWellModel &well_model,
    const Simulator &ebos_simulator,
    DeferredLogger &deferred_logger,
    WellState &well_state
) :
    well_model_{well_model},
    ebos_simulator_{ebos_simulator},
    deferred_logger_{deferred_logger},
    well_state_{well_state},
    schedule_{ebos_simulator.vanguard().schedule()},
    summary_state_{ebos_simulator.vanguard().summaryState()},
    report_step_idx_{ebos_simulator_.episodeIndex()},
    glo_{schedule_.glo(report_step_idx_)},
    debug{false}
{

}

/****************************************
 * Public methods in alphabetical order
 ****************************************/

template<typename TypeTag>
double
GasLiftGroupInfo<TypeTag>::
alqRate(const std::string &group_name)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    return group_rate.alq();
}

template<typename TypeTag>
void
GasLiftGroupInfo<TypeTag>::
initialize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);
    initializeGroupRatesRecursive_(group);
    std::vector<std::string> group_names;
    std::vector<double> group_efficiency;
    initializeWell2GroupMapRecursive_(
        group, group_names, group_efficiency, /*current efficiency=*/1.0);
}

template<typename TypeTag>
double
GasLiftGroupInfo<TypeTag>::
gasRate(const std::string &group_name)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    return group_rate.gasRate();
}

template<typename TypeTag>
std::optional<double>
GasLiftGroupInfo<TypeTag>::
gasTarget(const std::string &group_name)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    return group_rate.gasTarget();
}

template<typename TypeTag>
std::vector<std::pair<std::string,double>> &
GasLiftGroupInfo<TypeTag>::
getWellGroups(const std::string &well_name)
{
    assert(this->well_group_map_.count(well_name) == 1);
    return this->well_group_map_[well_name];
}

template<typename TypeTag>
std::optional<double>
GasLiftGroupInfo<TypeTag>::
maxAlq(const std::string &group_name)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    return group_rate.maxAlq();
}

template<typename TypeTag>
double
GasLiftGroupInfo<TypeTag>::
oilRate(const std::string &group_name)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    return group_rate.oilRate();
}

template<typename TypeTag>
std::optional<double>
GasLiftGroupInfo<TypeTag>::
oilTarget(const std::string &group_name)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    return group_rate.oilTarget();
}

template<typename TypeTag>
void
GasLiftGroupInfo<TypeTag>::
update(
    const std::string &group_name, double delta_oil, double delta_gas, double delta_alq)
{
    auto &group_rate = this->group_rate_map_.at(group_name);
    group_rate.update(delta_oil, delta_gas, delta_alq);
}


/****************************************
 * Private methods in alphabetical order
 ****************************************/


template<typename TypeTag>
bool
GasLiftGroupInfo<TypeTag>::
checkDoGasLiftOptimization_(const std::string &well_name)
{
    if (this->well_state_.gliftCheckAlqOscillation(well_name)) {
        displayDebugMessage_(
             "further optimization skipped due to oscillation in ALQ", well_name);
        return false;
    }
    if (this->optimize_only_thp_wells_) {
        //const auto &well = this->schedule_.getWell(well_name, this->report_step_idx_);
        WellInterfaceRawPtr well = this->well_model_.maybeGetWell(well_name);
        if (well) {
            const int well_index = well->indexOfWell();
            const Well::ProducerCMode& control_mode
                = this->well_state_.currentProductionControls()[well_index];
            if (control_mode != Well::ProducerCMode::THP ) {
                displayDebugMessage_("Not THP control. Skipping.", well_name);
                return false;
            }
        }
        else {
            // well_name is not present in the well_model's well container
            return false;
        }
    }
    if (!checkNewtonIterationIdxOk_(well_name)) {
        return false;
    }
    if (!this->glo_.has_well(well_name)) {
        displayDebugMessage_(
             "Gas Lift not activated: WLIFTOPT is probably missing", well_name);
        return false;
    }
    auto increment = this->glo_.gaslift_increment();
    // NOTE: According to the manual: LIFTOPT, item 1, :
    //   "Increment size for lift gas injection rate. Lift gas is
    //   allocated to individual wells in whole numbers of the increment
    //   size.  If gas lift optimization is no longer required, it can be
    //   turned off by entering a zero or negative number."
    if (increment <= 0) {
        if (this->debug) {
            const std::string msg = fmt::format(
                "Gas Lift switched off in LIFTOPT item 1 due to non-positive "
                "value: {}", increment);
                displayDebugMessage_(msg, well_name);
        }
        return false;
    }
    else {
        return true;
    }
}

template<typename TypeTag>
bool
GasLiftGroupInfo<TypeTag>::
checkNewtonIterationIdxOk_(const std::string &well_name)
{
    const int iteration_idx =
        this->ebos_simulator_.model().newtonMethod().numIterations();
    if (this->glo_.all_newton()) {
        const int nupcol = this->schedule_[this->report_step_idx_].nupcol();
        if (this->debug) {
            const std::string msg = fmt::format(
                "LIFTOPT item4 == YES, it = {}, nupcol = {} -->  GLIFT optimize = {}",
                iteration_idx,
                nupcol,
                ((iteration_idx <= nupcol) ? "TRUE" : "FALSE"));
            displayDebugMessage_(msg, well_name);
        }
        return iteration_idx <= nupcol;
    }
    else {
        if (this->debug) {
            const std::string msg = fmt::format(
                    "LIFTOPT item4 == NO, it = {} --> GLIFT optimize = {}",
                    iteration_idx, ((iteration_idx == 1) ? "TRUE" : "FALSE"));
            displayDebugMessage_(msg, well_name);
        }
        return iteration_idx == 1;
    }
}

template<typename TypeTag>
void
GasLiftGroupInfo<TypeTag>::
displayDebugMessage_(const std::string &msg)
{
    if (this->debug) {
        const std::string message = fmt::format(
             "  GLIFT (DEBUG) : Init group info : {}", msg);
        this->deferred_logger_.info(message);
    }
}

template<typename TypeTag>
void
GasLiftGroupInfo<TypeTag>::
displayDebugMessage_(const std::string &msg, const std::string &well_name)
{
    if (this->debug) {
        const std::string message = fmt::format(
             "  GLIFT (DEBUG) : Init group info : Well {} : {}",
             well_name, msg);
        this->deferred_logger_.info(message);
    }
}


template<typename TypeTag>
std::pair<double, double>
GasLiftGroupInfo<TypeTag>::
getProducerWellRates_(WellInterfaceRawPtr well)
{
    //const auto &well = this->schedule_.getWell(well_name, this->report_step_idx_);
    assert(well->isProducer());
    const int well_index = well->indexOfWell();
    const int np = this->well_state_.numPhases();
    const auto& pu = well->phaseUsage();
    auto oil_rate =
        -this->well_state_.wellRates()[np * well_index + pu.phase_pos[Oil]];
    auto gas_rate =
        -this->well_state_.wellRates()[np * well_index + pu.phase_pos[Gas]];
    return {oil_rate, gas_rate};
}

template<typename TypeTag>
void
GasLiftGroupInfo<TypeTag>::
initializeWell2GroupMapRecursive_(
    const Group &group,
    std::vector<std::string> &group_names,
    std::vector<double> &group_efficiency,
    double cur_efficiency)
{
    double gfac = group.getGroupEfficiencyFactor();
    cur_efficiency = gfac * cur_efficiency;
    for (auto &item : group_efficiency) {
        item *= gfac;
    }
    if (this->group_rate_map_.count(group.name()) == 1) {
        // extract the subset of groups that has limits or targets that can affect
        //   gas lift optimization.
        group_names.push_back(group.name());
        group_efficiency.push_back(gfac);
    }
    if (group.wellgroup()) {
        for (const std::string& well_name : group.wells()) {
            // TODO: can the same well be memember of two different groups
            //  (on the same recursion level) ?
            assert(this->well_group_map_.count(well_name) == 0);
            if (checkDoGasLiftOptimization_(well_name)) {
                const auto &well = this->schedule_.getWell(
                    well_name, this->report_step_idx_);
                double wfac = well.getEfficiencyFactor();
                auto [itr, success] = this->well_group_map_.insert(
                      {well_name, /*empty vector*/ {}});
                assert(success);
                auto &vec = itr->second;
                assert(group_names.size() == group_efficiency.size());
                auto iter2 = group_efficiency.begin();
                for (auto iter1 = group_names.begin();
                     iter1 != group_names.end(); ++iter1)
                {
                    double efficiency = (*iter2) * wfac;
                    vec.emplace_back(/*group_name=*/*iter1, efficiency);
                    ++iter2;
                }
            }
        }
    }
    else {
        for (const std::string& group_name : group.groups()) {
            if (!this->schedule_.back().groups.has(group_name))
                continue;
            const Group& sub_group = this->schedule_.getGroup(
                group_name, this->report_step_idx_);
            initializeWell2GroupMapRecursive_(
                sub_group, group_names, group_efficiency, cur_efficiency);
        }
    }
    if (this->group_rate_map_.count(group.name()) == 1) {
        group_names.pop_back();
        group_efficiency.pop_back();
    }
}

template<typename TypeTag>
std::tuple<double, double, double>
GasLiftGroupInfo<TypeTag>::
initializeGroupRatesRecursive_(const Group &group)
{
    double oil_rate = 0.0;
    double gas_rate = 0.0;
    double alq = 0.0;
    if (group.wellgroup()) {
        for (const std::string& well_name : group.wells()) {
            // NOTE: The well does not have to be present in the well_model's
            //   well container..
            WellInterfaceRawPtr well = this->well_model_.maybeGetWell(well_name);
            if (well) {
                if (well->isProducer()) {
                    auto [sw_oil_rate, sw_gas_rate] = getProducerWellRates_(well);
                    auto sw_alq = this->well_state_.getALQ(well_name);
                    //const auto &well_ecl =
                    //    this->schedule_.getWell(well_name, this->report_step_idx_);
                    double factor = well->wellEcl().getEfficiencyFactor();
                    oil_rate += (factor * sw_oil_rate);
                    gas_rate += (factor * sw_gas_rate);
                    alq += (factor * sw_alq);
                }
            }
        }
    }
    else {
        for (const std::string& group_name : group.groups()) {
            if (!this->schedule_.back().groups.has(group_name))
                continue;
            const Group& sub_group = this->schedule_.getGroup(
                group_name, this->report_step_idx_);
            auto [sg_oil_rate, sg_gas_rate, sg_alq]
                = initializeGroupRatesRecursive_(sub_group);
            const auto gefac = sub_group.getGroupEfficiencyFactor();
            oil_rate += (gefac * sg_oil_rate);
            gas_rate += (gefac * sg_gas_rate);
            alq += (gefac * sg_alq);
        }
    }
    std::optional<double> oil_target, gas_target, max_total_gas, max_alq;
    const auto controls = group.productionControls(this->summary_state_);
    if (group.has_control(Group::ProductionCMode::ORAT)) {
        oil_target = controls.oil_target;
    }
    if (group.has_control(Group::ProductionCMode::GRAT)) {
        gas_target = controls.gas_target;
    }
    if (this->glo_.has_group(group.name())) {
        const auto &gl_group = this->glo_.group(group.name());
        max_alq = gl_group.max_lift_gas();
        max_total_gas = gl_group.max_total_gas();
    }
    if (oil_target || gas_target || max_total_gas || max_alq) {
        this->group_rate_map_.try_emplace(group.name(),
            oil_rate, gas_rate, alq, oil_target, gas_target, max_total_gas, max_alq);
    }
    return std::make_tuple(oil_rate, gas_rate, alq);
}
