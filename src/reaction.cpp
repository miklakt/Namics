#include "reaction.h"

Reaction::Reaction(const Input* In_,vector<Segment*> Seg_, vector<State*> Sta_, string name_) {
	In=In_; name=name_;   Sta=Sta_; Seg=Seg_;
	KEYS.push_back("K"); 
	KEYS.push_back("pK");
	KEYS.push_back("equation");  
}
Reaction::~Reaction() = default;


void Reaction::PutParameter(string new_param) {
NAMICS_DBG("PutParameter in Reaction " + name << endl); KEYS.push_back(new_param); 
}

bool Reaction::CheckInput(int start) {
NAMICS_DBG("CheckInput in Reaction " + name << endl);	bool success=true;
	K=-1;
	pK=-100;
	Sto.clear();
	State_nr.clear();
	success= In->CheckParameters("reaction",name,start, KEYS, PARAMETERS);
	if (success) {
		if (GetValue("K").size()==0 && GetValue("pK").size()==0)  {
			cout <<" reaction " << name << " has no K nor pK value" << endl; success=false;
		} else {
			if (GetValue("K").size() ==0) {
				pK=ParseReal(GetValue("pK"),pK);
				if (pK==-100) {cout <<" reaction " << name << " no valid pK value found " << endl; success=false; }
				K=pow(10,-pK); 
			} else {
				K=ParseReal(GetValue("K"),pK);
				if (K<0) {
					cout <<" reaction " << name << " has not a positive value for 'K' " << endl; success=false;
				} else {
					pK=-log10(K); 
				}
			}
		}
		equation=GetValue("equation"); 
		if (equation.size()==0) {
			success=false; cout <<" reaction " << name << " has no equation specified" <<endl; 
		} else {
			string s=equation; 
			vector<string>sub_equal;
			In->split(s,'=',sub_equal);
			if (sub_equal.size()!=2) {
				cout <<" reaction : " << name << " equation : " << equation << "should have one '=' sign" << endl; 
				success=false; 
			} else {
				for (int k=0; k<2; k++) {
					vector<string>sub_plus;
					In->split(sub_equal[k],'+',sub_plus);
					int sub_l=sub_plus.size();
					for (int l=0; l<sub_l; l++) { 
						vector<int>open;
						vector<int>close;
						In->EvenBrackets(sub_plus[l], open, close);
						int length=open.size();
						if (length !=1) {
							cout <<" reaction : " << name << " equation " << equation << " has too many mon types in between '+' signs " << endl; 
							success=false; 
						} else  {
							string state_name=sub_plus[l].substr(open[0]+1,close[0]-open[0]-1);
							int num_states=In->StateList.size(); 
							bool found=false;
							for (int i=0; i<num_states; i++) {
								string s_name=In->StateList[i];
								if (state_name == s_name) {
									found = true; 
									State_nr.push_back(i); 
									Sta[i]->in_reaction=true;
									Seg_nr.push_back(Sta[i]->mon_nr);
									int state_length=Seg[Sta[i]->mon_nr]->state_name.size();
									for (int k=0; k<state_length; k++) 
										if (Seg[Sta[i]->mon_nr]->state_name[k]==state_name) State_in_seg_nr.push_back(k);
								}
							}
							if (!found) {cout << " reaction : " << name << " equation " << equation << " state " << state_name << " not found " << endl; success=false;  }
							
							int sto=ParseInt(sub_plus[l].substr(0,open[0]),0);
							if (sto<1) {
								if (sto==0) cout << " reaction : " << name << " equation : " << equation << " has a zero as stocheometry number " << endl;
								else
								cout << " reaction : " << name << " equation : " << equation << " has negative stocheometry number " << endl; 
								success=false; 
							}
							if (k==0) sto*=-1; //lhs terms get negative stocheometry numbers.
							Sto.push_back(sto);
						}
					}
				}
				
			}

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
			cout <<" reaction : " << name << " equation : " << equation << " not balanced for internal states of mon type: " << In->MonList[i] << endl; 
		}
	}
	Real charge=0;
	for (int j=0; j<LENGTH; j++) {
		charge += Sto[j]*Sta[State_nr[j]]->valence;
	}
	if (abs(charge)>1e-10) {
		success=false;
		cout <<" reaction : " << name << " equation : " << equation << " not electroneutral. Absolute value appears : " << abs(charge) << endl; 
	}
	return success;
}
 
string Reaction::GetValue(string parameter){
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return ""; 
}
 
void Reaction::push(string s, Real X) {
NAMICS_DBG("push (Real) in Reaction " + name << endl); Reals.push_back(s);
	Reals_value.push_back(X); 
}
void Reaction::push(string s, int X) {
NAMICS_DBG("push (int) in Reaction " + name << endl); ints.push_back(s);
	ints_value.push_back(X); 
}
void Reaction::push(string s, bool X) {
NAMICS_DBG("push (boool) in Reaction " + name << endl); bools.push_back(s);
	bools_value.push_back(X); 
}
void Reaction::push(string s, string X) {
NAMICS_DBG("push (string) in Reaction " + name << endl); strings.push_back(s);
	strings_value.push_back(X); 	
}
void Reaction::PushOutput() {
NAMICS_DBG("PushOutput in Reaction " + name << endl); strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();  
	push("equation",equation);
	push("pK",pK);
}

Real* Reaction::GetPointer(string s,int &SIZE) {
	(void)SIZE;
	(void)s;
NAMICS_DBG("GetPointer in Reaction " + name << endl);	return NULL;
}
int* Reaction::GetPointerInt(string s, int &SIZE) {
	(void)SIZE;
	(void)s;
NAMICS_DBG("GetPointerInt in Reaction " + name << endl);	return NULL;
}


int Reaction::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
NAMICS_DBG("GetValue (long)  in Reaction " + name << endl);	int i=0;
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

Real Reaction::ChemIntBulk(State* sta) {
	Real value=0;
	
	int mon_length=In->MonList.size();
	int state_length=In->StateList.size();
	for (int i=0; i<mon_length; i++) 
		if (Seg[i]->ns<2) {value+=sta->chi[i]*Seg[i]->phibulk;}
	for (int i=0; i<state_length; i++) {value+=sta->chi[mon_length+i]*Seg[Sta[i]->mon_nr]->state_phibulk[Sta[i]->state_nr];}
	return value;
}


Real Reaction::pKeff() {
	Real value=0;
	int length=Sto.size();
	for (int i=0; i<length; i++) {
		value +=Sto[i]*ChemIntBulk(Sta[i]);
	}
	return pK+value/log(10.0);
}

Real Reaction::Residual_value() { //only working when chi are not state dependent. chem is in log en not log10 but cal for pKeff should be done with ln..
	Real res_value=0;
	Real alphab=0;
	int length=Sto.size();
	for (int i=0; i<length; i++) {
		alphab=Seg[Seg_nr[i]]->state_alphabulk[State_in_seg_nr[i]];
		if (alphab>0) res_value-=Sto[i]*log10(alphab); else {
			
			cout <<"alphabulk =0 " <<  endl; 
			cout << "Seg " << Seg[Seg_nr[i]]->name << " State_in_seg_nr [i] " << State_in_seg_nr[i] << endl;
		}
	}
	return -1.0+res_value/pKeff();
}

bool Reaction::PutAlpha(Real alpha) {
	bool success=true;
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

	return success; 
}

bool Reaction::GuessAlpha() {
	bool success=true;
	Real alpha_f=-1;
	int length=Sto.size();
	Real k=pow(10,-pKeff());
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
			cout <<"Please fix alphabulk of a (charged) state of 'water'" << endl;
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
				k=k/pow(Seg[water]->state_alphabulk[State_in_seg_nr[i]],Sto[i]);
			} else  {
				if (state1<0) { 
					state1=State_in_seg_nr[i];
				} else {
					if (state2<0) state2=State_in_seg_nr[i];
					else {cout <<"more than 2 states found...." << endl; }
					if (Sto[i]>0) Seg[other]->ItState=State_in_seg_nr[i]; else
					cout <<"expected positive sto number for state2 " << endl;
				}	
			}
		}
		Seg[other]->state_alphabulk[state2]=k/(k+1);
		Seg[other]->state_alphabulk[state1]=1/(k+1); 



	}

 	return success;
}
