#include "state.h"

State::State(const Input* In_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_) {
	In=In_; name=name_;  Seg=Seg_;
}
State::~State() = default;

bool State::CheckInput(int start) {
NAMICS_DBG("CheckInput in State " + name << std::endl);	bool success=true;
	fixed=false;
	in_reaction =false;
	alphabulk =-1;
	valence=0;
	mon_nr=-1;
	state_id=-1;
	state_nr_of_copy =-1;
	seg_nr_of_copy =-1;
	int length=In->MonList.size();
	for (int k=0; k<length; k++) {
		if (Seg[k]->name == name) { std::cout << "name of state can not be the same as the name of any mon in the system" << std::endl;
		success=false;
		}
	}
	int length_StateList=In->StateList.size();
	for (int k=0; k<length_StateList; k++) {
		if (In->StateList[k] == name) state_id=k;
	}
	if (state_id<0) std::cout << "error: state_id is negative " << std::endl;
	const auto& parameters = In->Parameters("state", name, start);
	static const std::vector<std::string> keys = {"alphabulk", "valence", "mon"};
	chi.clear();
	success = true;
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		if (it.key().rfind("chi_", 0) == 0) {
			const std::string target = it.key().substr(4);
			if (ContainsValue(In->MonList, target) || ContainsValue(In->StateList, target)) {
				try {
					const Real chi_value = it.value().get<Real>();
					if (target == name && chi_value != 0) {
						std::cout <<" chi value for chi("<<name<<","<<target<<") value ignored: set to zero!" << std::endl;
					}
				} catch (const nlohmann::json::exception& error) {
					success = false;
					std::cout << "Invalid json type in state '" << name << "' for '" << it.key() << "': " << error.what() << std::endl;
				}
				continue;
			}
		}
		success = false;
		std::cout << "state property '" << it.key() << "' is unknown. Use a standard state key or chi_<mon/state name>." << std::endl;
	}
	if (success) {
		try {
		if (!parameters.contains("mon")) {
			std::cout << "Please specify the 'mon' for state " << name << std::endl;
			success=false;
		} else {
			const std::string mon_name = parameters.at("mon").get<std::string>();
			bool mon_found=false;
			for (int k=0; k<length; k++) {
				if (Seg[k]->name==mon_name) {
					mon_found=true;
					mon_nr=k;

				}
			}
			if (!mon_found) {
				std::cout <<"in state " << name << " mon name " << mon_name << " not found in system " << std::endl;
				success=false;
			}
			if (parameters.contains("alphabulk")) {
				fixed=true;
				alphabulk=parameters.at("alphabulk").get<Real>();
				if (alphabulk <0 || alphabulk > 1) {
					std::cout << "for state " << name << " value for alphabulk is out of range 0 ... 1 " << std::endl;
					success=false;
				}
			}
			valence = 0;
			if (parameters.contains("valence")) {
				valence=parameters.at("valence").get<Real>();
				if (valence <-10 || valence > 10) {
					std::cout << "for state " << name << " value for valence is out of range -10 ... 10 " << std::endl;
					success=false;
				}
			}
		}
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in state '" << name << "': " << error.what() << std::endl;
			success = false;
		}
		if (success) state_nr=Seg[mon_nr]->AddState(state_id,alphabulk,valence,fixed);
	}
	return success;
}

void State::PushOutput() {
NAMICS_DBG("PushOutput in State " + name << std::endl);
	OUTPUT = nlohmann::ordered_json::object();
	alphabulk=Seg[mon_nr]->state_alphabulk[state_nr];
	OUTPUT["alphabulk"] = alphabulk;
	OUTPUT["valence"] = valence;
	for (size_t i = 0; i < In->MonList.size() && i < chi.size(); i++) OUTPUT["chi_" + In->MonList[i]] = chi[i];
	for (size_t i = 0; i < In->StateList.size() && In->MonList.size() + i < chi.size(); i++) OUTPUT["chi_" + In->StateList[i]] = chi[In->MonList.size() + i];
	OUTPUT["phi"] = {{"profile", 0}};
}

std::span<Real> State::GetPointer(int profile) {
NAMICS_DBG("GetPointer in State " + name << std::endl);
	if (profile != 0) return {};
	const int M = Seg[mon_nr]->lat->M;
	return std::span<Real>(Seg[mon_nr]->phi_state).subspan(static_cast<size_t>(state_nr) * M, static_cast<size_t>(M));
}
