#include "reaction.h"

Reaction::Reaction(const Input* In_,std::span<const std::unique_ptr<Segment>> Seg_, std::span<const std::unique_ptr<State>> Sta_, std::string name_) {
	In=In_; name=name_;   Sta=Sta_; Seg=Seg_;
}

bool Reaction::CheckInput(int start) {
NAMICS_DBG("CheckInput in Reaction " + name << std::endl);	bool success=true;
	K=-1;
	pK=-100;
	Sto.clear();
	State_nr.clear();
	const auto& parameters = In->Parameters("reaction", name, start);
	static const std::vector<std::string> keys = {"K", "pK", "equation"};
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		success = false;
		std::cout << "reaction property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : keys) std::cout << item << std::endl;
	}
	if (success) {
		try {
		if (!parameters.contains("K") && !parameters.contains("pK"))  {
			std::cout <<" reaction " << name << " has no K nor pK value" << std::endl; success=false;
		} else {
			if (!parameters.contains("K")) {
				pK=parameters.at("pK").get<Real>();
				if (pK==-100) {std::cout <<" reaction " << name << " no valid pK value found " << std::endl; success=false; }
				K=std::pow(10,-pK); 
			} else {
				K=parameters.at("K").get<Real>();
				if (K<0) {
					std::cout <<" reaction " << name << " has not a positive value for 'K' " << std::endl; success=false;
				} else {
					pK=-log10(K); 
				}
			}
		}
		equation=parameters.value("equation", std::string{});
		if (equation.size()==0) {
			success=false; std::cout <<" reaction " << name << " has no equation specified" <<std::endl; 
		} else {
			std::string s=equation; 
			std::vector<std::string>sub_equal;
			In->split(s,'=',sub_equal);
			if (sub_equal.size()!=2) {
				std::cout <<" reaction : " << name << " equation : " << equation << "should have one '=' sign" << std::endl; 
				success=false; 
			} else {
				for (int k=0; k<2; k++) {
					std::vector<std::string>sub_plus;
					In->split(sub_equal[k],'+',sub_plus);
					int sub_l=sub_plus.size();
					for (int l=0; l<sub_l; l++) { 
						std::vector<int>open;
						std::vector<int>close;
						std::vector<char> stack;
						for (size_t i = 0; i < sub_plus[l].size(); ++i) {
							if (sub_plus[l][i] == '(') {
								stack.push_back(sub_plus[l][i]);
								open.push_back(static_cast<int>(i));
							} else if (sub_plus[l][i] == ')') {
								close.push_back(static_cast<int>(i));
								if (stack.empty()) break;
								stack.pop_back();
							}
						}
						int length=open.size();
						if (length !=1) {
							std::cout <<" reaction : " << name << " equation " << equation << " has too many mon types in between '+' signs " << std::endl; 
							success=false; 
						} else  {
							std::string state_name=sub_plus[l].substr(open[0]+1,close[0]-open[0]-1);
							int num_states=In->StateList.size(); 
							bool found=false;
							for (int i=0; i<num_states; i++) {
								std::string s_name=In->StateList[i];
								if (state_name == s_name) {
									found = true; 
									State_nr.push_back(i); 
									Seg_nr.push_back(Sta[i]->mon_nr);
									int state_length=Seg[Sta[i]->mon_nr]->state_name.size();
									for (int k=0; k<state_length; k++) 
										if (Seg[Sta[i]->mon_nr]->state_name[k]==state_name) State_in_seg_nr.push_back(k);
								}
							}
							if (!found) {std::cout << " reaction : " << name << " equation " << equation << " state " << state_name << " not found " << std::endl; success=false;  }
							
							int sto=ParseInt(sub_plus[l].substr(0,open[0]),0);
							if (sto<1) {
								if (sto==0) std::cout << " reaction : " << name << " equation : " << equation << " has a zero as stocheometry number " << std::endl;
								else
								std::cout << " reaction : " << name << " equation : " << equation << " has negative stocheometry number " << std::endl; 
								success=false; 
							}
							if (k==0) sto*=-1; //lhs terms get negative stocheometry numbers.
							Sto.push_back(sto);
						}
					}
				}
				
			}

		} 
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in reaction '" << name << "': " << error.what() << std::endl;
			success = false;
		}
	}
	int length=In->MonList.size();
	int LENGTH=Sto.size();  
	for (int i=0; i<length; i++) {
		int sum=0;
		
		for (int j=0; j<LENGTH; j++) {
			if (Seg_nr[j]==i) sum+=Sto[j]; 
		}
		if (sum !=0 ) {
			success=false;
			std::cout <<" reaction : " << name << " equation : " << equation << " not balanced for internal states of mon type: " << In->MonList[i] << std::endl; 
		}
	}
	Real charge=0;
	for (int j=0; j<LENGTH; j++) {
		charge += Sto[j]*Sta[State_nr[j]]->valence;
	}
	if (std::abs(charge)>1e-10) {
		success=false;
		std::cout <<" reaction : " << name << " equation : " << equation << " not electroneutral. Absolute value appears : " << std::abs(charge) << std::endl; 
	}
	return success;
}

void Reaction::PushOutput() {
NAMICS_DBG("PushOutput in Reaction " + name << std::endl);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["equation"] = equation;
	OUTPUT["pK"] = pK;
}

Real Reaction::pKeff() {
	Real value=0;
	int length=Sto.size();
	for (int i=0; i<length; i++) {
		const State& sta = *Sta[i];
		Real chem_bulk = 0;
		int mon_length = In->MonList.size();
		int state_length = In->StateList.size();
		for (int j=0; j<mon_length; j++)
			if (Seg[j]->ns<2) chem_bulk += sta.chi[j]*Seg[j]->phibulk;
		for (int j=0; j<state_length; j++)
			chem_bulk += sta.chi[mon_length+j]*Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
		value += Sto[i]*chem_bulk;
	}
	return pK+value/std::log(10.0);
}

Real Reaction::Residual_value() { //only working when chi are not state dependent. chem is in log en not log10 but cal for pKeff should be done with ln..
	Real res_value=0;
	Real alphab=0;
	int length=Sto.size();
	for (int i=0; i<length; i++) {
		alphab=Seg[Seg_nr[i]]->state_alphabulk[State_in_seg_nr[i]];
		if (alphab>0) res_value-=Sto[i]*log10(alphab); else {
			
			std::cout <<"alphabulk =0 " <<  std::endl; 
			std::cout << "Seg " << Seg[Seg_nr[i]]->name << " State_in_seg_nr [i] " << State_in_seg_nr[i] << std::endl;
		}
	}
	return -1.0+res_value/pKeff();
}

void Reaction::PutAlpha(Real alpha) {
	int water=-1;
	int other=-1;
	int length=Sto.size();
	for (int i=0; i<length; i++) {
		if (!Seg[Seg_nr[i]]->state_change[State_in_seg_nr[i]]) water=Seg_nr[i]; 
	}
	for (int i=0; i<length; i++) {
		if (Seg_nr[i]!=water) other =Seg_nr[i];
	}
	if (other>-1) Seg[other]->PutAlpha(alpha); else Seg[water]->PutAlpha(alpha);
}

bool Reaction::GuessAlpha() {
	bool success=true;
	Real alpha_f=-1;
	int length=Sto.size();
	Real k=std::pow(10,-pKeff());
	int water=-1;
	int other=-1;
	int ns; 
	Real sum_alpha=0;	
	int state1=-1,state2=-1;
		

	if (length==3) {
		for (int i=0; i<length; i++) {
			if (!Seg[Seg_nr[i]]->state_change[State_in_seg_nr[i]]) {
				water =Seg_nr[i];
				alpha_f = Seg[Seg_nr[i]]->state_alphabulk[State_in_seg_nr[i]];
			}
		}

		ns=Seg[water]->ns; 
		if (ns!=3) {
			std::cout <<"Please fix alphabulk of a (charged) state of 'water'" << std::endl;
			success=false;
		}	
	
		for (int i=0; i<ns; i++) {
			if (Seg[water]->state_change[i]) {
				if(Seg[water]->state_valence[i]==0) Seg[water]->state_alphabulk[i]=0; else {
					Seg[water]->state_alphabulk[i]=k/alpha_f;
					Seg[water]->ItState=i; 
				}
			}
			sum_alpha+=Seg[water]->state_alphabulk[i]; 
		}
		for (int i=0; i<ns; i++) {
			if (Seg[water]->state_change[i]) {
				if(Seg[water]->state_valence[i]==0) {
					Seg[water]->state_alphabulk[i]=1-sum_alpha; 
				}
			}
		}


	} else {
		for (int i=0; i<length; i++) {
			if (!Seg[Seg_nr[i]]->state_change[State_in_seg_nr[i]]) water =Seg_nr[i]; 
		}
		
		for (int i=0; i<length; i++) {
			if (Seg_nr[i]!=water) other =Seg_nr[i]; 
		}


		for (int i=0; i<length; i++) {
			if (Seg_nr[i]==water) {
				k=k/std::pow(Seg[water]->state_alphabulk[State_in_seg_nr[i]],Sto[i]);
			} else  {
				if (state1<0) { 
					state1=State_in_seg_nr[i];
				} else {
					if (state2<0) state2=State_in_seg_nr[i];
					else {std::cout <<"more than 2 states found...." << std::endl; }
					if (Sto[i]>0) Seg[other]->ItState=State_in_seg_nr[i]; else
					std::cout <<"expected positive sto number for state2 " << std::endl;
				}	
			}
		}
		Seg[other]->state_alphabulk[state2]=k/(k+1);
		Seg[other]->state_alphabulk[state1]=1/(k+1); 



	}

 	return success;
}
