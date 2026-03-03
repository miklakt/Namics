
#include "output.h"

Output::Output(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_,vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_,System* Sys_,Solve_scf* New_,string name_,int outnr,int N_out) {
NAMICS_DBG("constructor in Output "<< endl);	In=In_; Lat = Lat_; Seg=Seg_; Sta=Sta_; Rea=Rea_; Mol=Mol_; Sys=Sys_; name=name_; n_output=N_out; output_nr=outnr;  New=New_;
	writer = io::SharedWriter();
	//KEYS.push_back("write_output");
	lat=Lat;
	KEYS.push_back("write_bounds");
	KEYS.push_back("append");
	KEYS.push_back("use_output_folder");
	KEYS.push_back("write");
	KEYS.push_back("clear");
	KEYS.push_back("DOS");
	KEYS.push_back("header_separator");
	KEYS.push_back("filename");
	input_error=false;
	bin_folder = "bin"; // folder in Namics where the binary is located
	use_output_folder = true; // LINUX ONLY, when you remove this, add it as a default to its CheckInputs part.
	n_starts = In->GetNumStarts();
	first=0;

}
Output::~Output() {
NAMICS_DBG("destructor in output " << endl);}
void Output::PutParameter(string new_param) {
NAMICS_DBG("PutParameter in Output " << endl); KEYS.push_back(new_param);
}

bool Output::Load() {
NAMICS_DBG("Load in output " << endl);	bool success=true;
	int molnr=0;
	success= In->LoadItems(name, OUT_key, OUT_name, OUT_prop);
	if (success) {
		int length=OUT_key.size();

		for (int i=0; i<length; i++) {
			if (OUT_key[i]=="mol"){
				vector<string> sub;
				In->split(OUT_prop[i],'*',sub);


				if (!(sub[0]==OUT_prop[i])){
					bool wildmon;
					if (sub[0] =="") wildmon=false; else wildmon=true;

					int k=0; int mollength=In->MolList.size();
					while (k<mollength) {
						if (In->MolList[k]==OUT_name[i]) molnr=k;
						k++;
					}
					if (wildmon) {
						int monlength=Mol[molnr]->MolMonList.size();
						for (int j=0; j<monlength; j++) {
							OUT_key.push_back(OUT_key[i]);
							OUT_name.push_back(OUT_name[i]);
							string s=sub[0];
							s=s.append(Seg[Mol[molnr]->MolMonList[j]]->name);
							s=s.append(sub[1]);
							OUT_prop.push_back(s);
						}
					} else {
						cout << "Alias-based wildcard output is not supported in this minimal build." << endl;
						return false;
					}
					OUT_key.erase(OUT_key.begin()+i);
					OUT_name.erase(OUT_name.begin()+i);
					OUT_prop.erase(OUT_prop.begin()+i);

				}

			}
		}

	}
	return success;
}

bool Output::CheckInput(int start_) {
NAMICS_DBG("CheckInput in output " << endl);	start=start_;
	if (name != "kal" && name != "pro") {
		write = false;
		cout << "Output type '" << name << "' is disabled. Only 'kal' and 'pro' are supported." << endl;
		return true;
	}

	bool success=true;
	success=In->CheckParameters("output",name,start, KEYS, PARAMETERS);
	if (success) {
		DOS=false;
		if (GetValue("DOS").size()>0) {
			DOS=In->Get_bool(GetValue("DOS"),DOS);
		}

		if (GetValue("append").size()>0) {
			append=In->Get_bool(GetValue("append"),append);

			if (name=="pro") {
					if (append) cout << "Warning: for output of type 'pro', the append is set to 'false'." << endl;
			}
			if (first==0) first=start;
		} else {
			if (name=="kal") append=false;
			if (name=="pro") append=false;
		}

		write_bounds = In->Get_bool(GetValue("write_bounds"),false);
		write  = In->Get_bool(GetValue("write"),true);

		if (GetValue("header_separator").size()>0) {
			sep=GetValue("header_separator");
			if (sep=="classic") {
				sep = ":";
			} else {
				if (sep.length()>1) {
					cout <<"For output entry 'header_separator', expected to find the keyword 'classic' (meaning ':') or a single character." << endl;
					cout <<"header_separator set to the default value '_'" << endl;
					sep="_";
				}
			}
		} else {
			sep = "_";
		}

		if (GetValue("use_output_folder").size()>0) {
			use_output_folder = In->Get_bool(GetValue("use_output_folder"),use_output_folder);
		} // default is set in the constructor

		if (success) {
			if (!Load()) {
				cout <<"Error in Load() in output" << endl;
				success=false;
			}
		}
		write_option="no_error";
		if (GetValue("write_output").size()>0) {
			vector<string> option_list;
			option_list.push_back("always");
			option_list.push_back("no_error");
			option_list.push_back("never");
			if (!In->Get_string(GetValue("write_output"),write_option,option_list,"In output: 'write_output' not recognised. Use 'always', 'never', or 'no_error'. The last value is default.")){
				cout <<"continue with write_output : no_error" << endl;
			}
		}
	} else cout <<"Error in CheckParameters in output" << endl;
	return success;
}

string Output::GetValue(string parameter) {
NAMICS_DBG("GetValue in output " << endl); auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

int* Output::GetPointerInt(string key, string name, string prop, int &Size) {
NAMICS_DBG("GetPointerInt in output " << endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int listlength;
	int choice;
	int i,j;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if (key=="lat") choice=4;
	if (key=="output") choice = 5;

	switch(choice) {
		case 1:
			listlength=Sys->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Sys->strings[j]) return Sys->GetPointerInt(Sys->strings_value[j],Size);
				j++;
			}
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In->MolList[i]) {
					listlength= Mol[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointerInt(Mol[i]->strings_value[j],Size);
						j++;
					}
				}
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In->MonList[i]) {
					listlength= Seg[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointerInt(Seg[i]->strings_value[j],Size);
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			listlength= lat->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==lat->strings[j]) return lat->GetPointerInt(lat->strings_value[j],Size);
				j++;
			}
			break;
		case 5:
			listlength=PointerVectorInt.size();
			j=0;
			while (j<listlength) {
				if (prop==strings[j]){ Size=SizeVectorInt[j];  return PointerVectorInt[j];}
				j++;
			}
			break;
		default:
			cout << "Program error: in Output, GetPointerInt reaches default...." << endl;
	}
	return NULL;
}
Real* Output::GetPointer(string key, string name, string prop, int &Size) {
NAMICS_DBG("GetPointer in output " << endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int listlength;
	int choice;
	int i,j;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if (key=="lat") choice=4;
	if (key=="output") choice = 5;

	switch(choice) {
		case 1:
			listlength=Sys->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Sys->strings[j]) return Sys->GetPointer(Sys->strings_value[j],Size);
				j++;
			}
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In->MolList[i]) {
					listlength= Mol[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointer(Mol[i]->strings_value[j],Size);
						j++;
					}
				}
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In->MonList[i]) {
					listlength= Seg[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointer(Seg[i]->strings_value[j],Size);
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			listlength= lat->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==lat->strings[j]) return lat->GetPointer(lat->strings_value[j],Size);
				j++;
			}
			break;
		case 5:
			listlength=PointerVectorReal.size();
			j=0;
			while (j<listlength) {
				if (prop==strings[j]){ Size=SizeVectorReal[j];  return PointerVectorReal[j];}
				j++;
			}
			break;
		default:
			cout << "Program error: in Output, GetPointer reaches default...." << endl;
	}
	return NULL;
}
int Output::GetValue(string key, string name, string prop, int &int_result, Real &Real_result, string &string_result) {
NAMICS_DBG("GetValue (long) in output " << endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int choice=0;
	int i;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if  (key=="newton") choice=4;
	if  (key=="lat") choice=5;
	if  (key=="output") choice=6;
	switch(choice) {
		case 1:
			return Sys->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In->MolList[i]) return Mol[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In->MonList[i]) return Seg[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 4:
			return New->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 5:
			return lat->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 6:
			return GetValue(prop,name,int_result,Real_result,string_result);
			break;
		default:
			cout << "Program error: in Output, GetValue reaches default...." << endl;
	}
	return 0;
}

void Output::WriteOutput(int subl) {
NAMICS_DBG("WriteOutput in output " + name << endl);	lat->subl=subl;
	if (!write) return;
	if (name != "kal" && name != "pro") {
		return;
	}
	int Size=0;
	string filename;
	vector<string> sub;

	string infilename = In->name;
	if (GetValue("filename").size()>0) infilename=GetValue("filename");
	In->split(infilename,'.',sub);
	string key;

	if (use_output_folder == true) {

		int occurrences = 0;
		string::size_type start = 0;
		string slash = "/";

		// Check if we have any slashes in the filename (are we in inputs, or higher up?)
		while ((start = sub[0].find(slash, start)) != string::npos) {
    	++occurrences;
    	start += slash.length();
		}

		// If we're not in the inputs folder, discard the path and take only the filename
		if  (occurrences != 0) {
			size_t found = sub[0].find_last_of("/\\");
			sub[0] = sub[0].substr(found+1);
		}
	}

    string numc = to_string(subl);
    string numcc = to_string(start);

    if (name=="kal") filename=sub[0].append(".").append(name); else {
		if (n_starts==1 && subl < 1) filename=sub[0].append(".").append(name);
		if (n_starts==1 && subl >0) filename = sub[0].append("_").append(numc).append(".").append(name);
		if (n_starts>1  && subl < 1) filename = sub[0].append("_").append(numcc).append(".").append(name);
		if (n_starts>1 && subl >0)  filename=sub[0].append("_").append(numcc).append("_").append(numc).append(".").append(name);
	}

	filename = In->output_info.getOutputPath() + filename;
	if (name=="pro") {
		vector<Real*> pointer;
		FILE *fp;
		fp=writer->OpenRaw(filename.c_str(),"w");
		int length=OUT_key.size();
		switch(lat->gradients) {
			case 1:
				writer->Writef(fp,"x\t");
				break;
			case 2:
				writer->Writef(fp,"x\ty\t");
				break;
			case 3:
				writer->Writef(fp,"x\ty\tz\t");
				break;
			default:
				break;
		}
		for (int i=0; i<length; i++) {
			string s=OUT_key[i];
			s.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);
			Real*  X = GetPointer(OUT_key[i],OUT_name[i],OUT_prop[i],Size);
			if (X!=NULL) {
				pointer.push_back(X);
				key = OUT_key[i];
				s = key.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);
				if (i<length-1) writer->Writef(fp,"%s\t",s.c_str()); else writer->Writef(fp,"%s",s.c_str());
			} else {cout << " Error for 'pro' output. It is only possible to output quantities known to be a 'profile'. That is why output quantity " + s + " is rejected. " << endl;}
		}
		if (DOS) writer->Writef(fp,"\r\n"); else writer->Writef(fp,"\n");
		Lat-> PutProfiles(fp,pointer,write_bounds,DOS);

		writer->Close(fp);
	}

	if (name=="kal") {
		if (start>first || subl>0) append=true;
		FILE *fp;
		if (!(writer->Exists(filename) && append)) {
			fp=writer->OpenRaw(filename.c_str(),"w");
			int length = OUT_key.size();
			for (int i=0; i<length; i++) {
				key=OUT_key[i];
				string s=key.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);
				if (i<length-1) writer->Writef(fp,"%s\t",s.c_str()); else writer->Writef(fp,"%s",s.c_str());
			}
			if (DOS) writer->Writef(fp,"\r\n"); else writer->Writef(fp,"\n");
		} else fp=writer->OpenRaw(filename.c_str(),"a");

		if (fp == NULL) {
			cerr << "Error trying to open " << filename.c_str() << endl;
			perror("Error");
		}


		int length = OUT_key.size();
		for (int i=0; i<length; i++) {
			int int_result=0;
			int result_nr=0;
			Real Real_result=0;
			string string_result;
			vector<string> sub;
			In-> split(OUT_prop[i],'(',sub);
			result_nr= GetValue(OUT_key[i],OUT_name[i],sub[0],int_result,Real_result,string_result);
			if (result_nr==0) {if (i<length-1) writer->Writef(fp,"NiN\t"); else writer->Writef(fp,"NiN");}
			if (result_nr==1) {if (i<length-1) writer->Writef(fp,"%i\t",int_result); else writer->Writef(fp,"%i",int_result);}
#ifdef LongReal
			if (result_nr==2) {if (i<length-1) writer->Writef(fp,"%.16Le\t",Real_result); else  writer->Writef(fp,"%.16Le",Real_result);}
#else
			if (result_nr==2) {if (i<length-1) writer->Writef(fp,"%.16e\t",Real_result); else  writer->Writef(fp,"%.16e",Real_result);}
#endif
			if (result_nr==3) {
				if (sub[0]==OUT_prop[i]) {
					if (i<length-1) writer->Writef(fp,"%s\t",string_result.c_str()); else writer->Writef(fp,"%s",string_result.c_str());
				} else {
					Real* X=GetPointer(OUT_key[i],OUT_name[i],sub[0],Size);
#ifdef LongReal
					if (i<length-1) writer->Writef(fp,"%.16Le\t",lat->GetValue(X,sub[1])); else writer->Writef(fp,"%.16Le",lat->GetValue(X,sub[1]));
#else
					if (i<length-1) writer->Writef(fp,"%.16e\t",lat->GetValue(X,sub[1])); else writer->Writef(fp,"%.16e",lat->GetValue(X,sub[1]));
#endif
				}
			}
		}
		if (DOS) writer->Writef(fp,"\r\n"); else writer->Writef(fp,"\n");
		writer->Close(fp);
	}
}

int Output::GetValue(string prop, string mod, int& int_result, Real& Real_result, string& string_result) {
	(void)mod;
  int i = 0;
  int_result=0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  Real_result=0.0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  string_result="false";
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}

void Output::push(string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Output::push(string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Output::push(string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}

void Output::push(string s, string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
