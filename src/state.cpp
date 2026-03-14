#include "state.h"

State::State(const Input* In_,vector<Segment*> Seg_, string name_) {
	In=In_; name=name_;  Seg=Seg_;
	KEYS.push_back("alphabulk");
	KEYS.push_back("valence");
	KEYS.push_back("mon");
}
State::~State() = default;


void State::PutParameter(string new_param) {
NAMICS_DBG("PutParameter in State " + name << endl); KEYS.push_back(new_param);
}

bool State::CheckInput(int start) {
NAMICS_DBG("CheckInput in State " + name << endl);	bool success=true;
	fixed=false;
	in_reaction =false;
	alphabulk =-1;
	valence=0;
	state_id=-1;
	state_nr_of_copy =-1;
	seg_nr_of_copy =-1;
	success= In->CheckParameters("state",name,start, KEYS, PARAMETERS);
	int length=In->MonList.size();
	for (int k=0; k<length; k++) {
		if (Seg[k]->name == name) { cout << "name of state can not be the same as the name of any mon in the system" << endl;
		success=false;
		}
	}
	int length_StateList=In->StateList.size();
	for (int k=0; k<length_StateList; k++) {
		if (In->StateList[k] == name) state_id=k;
	}
	if (state_id<0) cout << "error: state_id is negative " << endl;
	if (success) {
		if (GetValue("mon").size()==0) {
			cout << "Please specify the 'mon' for state " << name << endl;
			success=false;
		} else {
			string mon_name;
			mon_name=GetValue("mon");
			bool mon_found=false;
			for (int k=0; k<length; k++) {
				if (Seg[k]->name==mon_name) {
					mon_found=true;
					mon_nr=k;

				}
			}
			if (!mon_found) {
				cout <<"in state " << name << " mon name " << mon_name << " not found in system " << endl;
				success=false;
			}
			if (GetValue("alphabulk").size()>0) {
				fixed=true;
				alphabulk=ParseReal(GetValue("alphabulk"),alphabulk);
				if (alphabulk <0 || alphabulk > 1) {
					cout << "for state " << name << " value for alphabulk is out of range 0 ... 1 " << endl;
					success=false;
				}
			}
			valence = 0;
			if (GetValue("valence").size()>0) {
				const string valence_value = GetValue("valence");
				valence=ParseReal(valence_value,valence);
				if (valence <-10 || valence > 10) {
					cout << "for state " << name << " value for valence " << valence_value << " is out of range -10 ... 10 " << endl;
					success=false;
				}
			}
		}
		state_nr=Seg[mon_nr]->AddState(state_id,alphabulk,valence,fixed);
	}
	length=chi_name.size();
	Real Chi;
	for (int i=0; i<length; i++) {
		Chi=-999;
		const string chi_value = GetValue("chi_"+chi_name[i]);
		if (chi_value.size()>0) {
			Chi=ParseReal(chi_value,Chi);
			if (Chi==-999) {success=false; cout <<" chi value: chi("<<name<<","<<chi_name[i]<<") = "<<chi_value << "not valid." << endl; }
			if (name==chi_name[i] && Chi!=0) {if (Chi!=-999) cout <<" chi value for chi("<<name<<","<<chi_name[i]<<") = "<<chi_value << "value ignored: set to zero!" << endl; Chi=0;}

		}
		chi[i]=Chi;
	}
	return success;
}

string State::GetValue(string parameter){
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

void State::PutChiKEY(string new_name) {
NAMICS_DBG("PutChiKey " + name << endl); KEYS.push_back("chi_" + new_name);
	chi_name.push_back(new_name);
	chi.push_back(-999);
}

void State::push(string s, Real X) {
NAMICS_DBG("push (Real) in State " + name << endl); Reals.push_back(s);
	Reals_value.push_back(X);
}
void State::push(string s, int X) {
NAMICS_DBG("push (int) in State " + name << endl); ints.push_back(s);
	ints_value.push_back(X);
}
void State::push(string s, bool X) {
NAMICS_DBG("push (boool) in State " + name << endl); bools.push_back(s);
	bools_value.push_back(X);
}
void State::push(string s, string X) {
NAMICS_DBG("push (string) in State " + name << endl); strings.push_back(s);
	strings_value.push_back(X);
}
void State::PushOutput() {
NAMICS_DBG("PushOutput in State " + name << endl); strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	alphabulk=Seg[mon_nr]->state_alphabulk[state_nr];
	push("alphabulk",alphabulk);
	push("valence",valence);
	int length=chi_name.size();
	for (int i=0; i<length; i++) push("chi_"+chi_name[i],chi[i]);
}

Real* State::GetPointer(string s,int &SIZE) {
	(void)SIZE;
	(void)s;
NAMICS_DBG("GetPointer in State " + name << endl);	return NULL;
}

int* State::GetPointerInt(string s,int &SIZE) {
	(void)SIZE;
	(void)s;
NAMICS_DBG("GetPointerInt in State " + name << endl);	return NULL;
}


int State::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
NAMICS_DBG("GetValue (long)  in State " + name << endl);	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true"; else string_result="false";
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}
		i++;
	}
	return 0;
}
