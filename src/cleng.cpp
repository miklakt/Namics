#include "cleng.h"
#include <random>

Cleng::Cleng(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_,  string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sys{Sys_},
      New{New_}
{
	if (debug) cout << "Cleng initialized" << endl;

    KEYS.push_back("MCS");
  	KEYS.push_back("save_interval");
  	KEYS.push_back("save_filename");
  	KEYS.push_back("seed");
}

Cleng::~Cleng() {
}

bool Cleng::CheckInput(int start) {
    if (debug) cout << "CheckInput in Cleng" << endl;
    bool success = true;
    
    success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
    if (success) {
        vector<string> options;
        if (GetValue("MCS").size() > 0) {
            success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of timesteps should be between 1 and 10000");
        }
        if (debug) cout << "MCS is " << MCS << endl;
        
        if (GetValue("save_interval").size() > 0) {
            success = In[0]->Get_int(GetValue("save_interval"), save_interval,1,MCS,"The save interval nr should be between 1 and 100");
        }
        if (debug) cout << "Save_interval " << save_interval << endl;
        
        if (Sys[0]->SysClampList.size() <1) {
            cout <<"Cleng needs to have clamped molecules in the system" << endl; success=false;
        }
        else {
            clamp_seg=Sys[0]->SysClampList[0]; 
            if (Sys[0]->SysClampList.size()>1) {
                success=false; cout <<"Currently the clamping is limited to one molecule per system. " << endl; 
            }
        }
        
        if (success) {
            n_boxes = Seg[clamp_seg]->n_box;
            sub_box_size=Seg[clamp_seg]->mx;
        }
        clp_mol=-1;
        int length = In[0]->MolList.size();
        for (int i=0; i<length; i++) if (Mol[i]->freedom =="clamped") clp_mol=i; 
    }
    if (success) {
        n_out = In[0]->OutputList.size();
        if (n_out == 0) cout << "Warning: no output defined!" << endl;
        
        for (int i =  0; i < n_out; i++) {
            Out.push_back(new Output(In, Lat, Seg, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
            if (!Out[i]->CheckInput(start)) {
                cout << "input_error in output " << endl;
                success=false;
            }
        }
        MonteCarlo();
    }
    return success;
}

bool Cleng::CP(transfer tofrom) {
    if (debug) cout << "CP in Cleng" << endl;
    int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int M = Lat[0]->M;
	bool success=true;
	int i,j;
	int length; 
	bool found;
	switch(tofrom) {
        case to_cleng:
			for (i=0; i<n_boxes; i++) {

			    PX=Seg[clamp_seg]->px1[i]; if (PX>MX) {PX-=MX; Sx.push_back(1);} else {Sx.push_back(0);}
				PY=Seg[clamp_seg]->py1[i]; if (PY>MY) {PY-=MY; Sy.push_back(1);} else {Sy.push_back(0);}
				PZ=Seg[clamp_seg]->pz1[i]; if (PZ>MZ) {PZ-=MZ; Sz.push_back(1);} else {Sz.push_back(0);}

				length=X.size();
				found=false; j=0;

				while (!found && j<length) {
					if (X[j]==PX && Y[j]==PY && Z[j]==PZ) {
					    found = true;
					    P.push_back(j);
					} else j++;
				}

				if (!found) {
					X.push_back(PX);
					Y.push_back(PY);
					Z.push_back(PZ);
					P.push_back(X.size()-1);
				}

				PX=Seg[clamp_seg]->px2[i]; if (PX>MX) {PX-=MX; Sx.push_back(1);} else {Sx.push_back(0);}
				PY=Seg[clamp_seg]->py2[i]; if (PY>MY) {PY-=MY; Sy.push_back(1);} else {Sy.push_back(0);}
				PZ=Seg[clamp_seg]->pz2[i]; if (PZ>MZ) {PZ-=MZ; Sz.push_back(1);} else {Sz.push_back(0);}
				length=X.size();

				found=false; j=0;
				while (!found && j<length) {
					if (X[j]==PX && Y[j]==PY && Z[j]==PZ) {
					    found = true;
					    P.push_back(j);
					}
					else j++;
				}
				if (!found) {
					X.push_back(PX);
					Y.push_back(PY);
					Z.push_back(PZ);
					P.push_back(X.size()-1);
				}
            }
			break;

        case to_segment:
			Zero(Seg[clamp_seg]->H_MASK,M);
			for (i=0; i < n_boxes; i++) {
				Seg[clamp_seg]->px1[i] = X[P[2*i]] + Sx[2*i]*MX;
                Seg[clamp_seg]->py1[i] = Y[P[2*i]] + Sy[2*i]*MY;
                Seg[clamp_seg]->pz1[i] = Z[P[2*i]] + Sz[2*i]*MZ;

				Seg[clamp_seg]->px2[i] = X[P[2*i+1]] + Sx[2*i+1]*MX;
				Seg[clamp_seg]->py2[i] = Y[P[2*i+1]] + Sy[2*i+1]*MY;
				Seg[clamp_seg]->pz2[i] = Z[P[2*i+1]] + Sz[2*i+1]*MZ;

				Seg[clamp_seg]->bx[i] = (Seg[clamp_seg]->px2[i] + Seg[clamp_seg]->px1[i] - sub_box_size) / 2;
				Seg[clamp_seg]->by[i] = (Seg[clamp_seg]->py2[i] + Seg[clamp_seg]->py1[i] - sub_box_size) / 2;
				Seg[clamp_seg]->bz[i] = (Seg[clamp_seg]->pz2[i] + Seg[clamp_seg]->pz1[i] - sub_box_size) / 2;

				if (Seg[clamp_seg]->bx[i] < 1) {Seg[clamp_seg]->bx[i] += MX; Seg[clamp_seg]->px1[i] +=MX; Seg[clamp_seg]->px2[i] +=MX;}
				if (Seg[clamp_seg]->by[i] < 1) {Seg[clamp_seg]->by[i] += MY; Seg[clamp_seg]->py1[i] +=MY; Seg[clamp_seg]->py2[i] +=MY;}
				if (Seg[clamp_seg]->bz[i] < 1) {Seg[clamp_seg]->bz[i] += MZ; Seg[clamp_seg]->pz1[i] +=MZ; Seg[clamp_seg]->pz2[i] +=MZ;}

				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px1[i]-1)%MX+1)*JX + ((Seg[clamp_seg]->py1[i]-1)%MY+1)*JY + (Seg[clamp_seg]->pz1[i]-1)%MZ+1]=1;
				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px2[i]-1)%MX+1)*JX + ((Seg[clamp_seg]->py2[i]-1)%MY+1)*JY + (Seg[clamp_seg]->pz2[i]-1)%MZ+1]=1;
            }
		break;

        default:
			success=false; 
			cout <<"error in transfer" << endl;
		break;
	}
    return success; 
}

void Cleng::WriteOutput(int subloop){
    if (debug) cout << "WriteOutput in Cleng" << endl;
      PushOutput();
      Sys[0]->PushOutput(); // needs to be after pushing output for seg.
      Lat[0]->PushOutput();
      New[0]->PushOutput();
      int length = In[0]->MonList.size();
      for (int i = 0; i < length; i++) Seg[i]->PushOutput();
      length = In[0]->MolList.size();
      for (int i = 0; i < length; i++) {
          int length_al = Mol[i]->MolAlList.size();
          for (int k = 0; k < length_al; k++) {
              Mol[i]->Al[k]->PushOutput();
          } 
          Mol[i]->PushOutput();
      }
      // length = In[0]->AliasList.size();
      // for (int i=0; i<length; i++) Al[i]->PushOutput();

      for (int i = 0; i < n_out; i++) {
          Out[i]->WriteOutput(subloop);
      }
}


int Cleng::GetIntRandomValue(int min_value, int max_value) {
    if (debug) cout << "Int GetRandomValue in Cleng" << endl;
    random_device rd;
    default_random_engine gen(rd());
    uniform_int_distribution<> dist(min_value, max_value);
    return dist(gen);
}

Real Cleng::GetRealRandomValue(int min_value, int max_value) {
    if (debug) cout << "Real GetRandomValue in Cleng" << endl;
    random_device rd;
    default_random_engine gen(rd());
    uniform_real_distribution<> dist(0, 1);
    return dist(gen);
}

bool Cleng::MakeShift() {
    if (debug) cout << "MakeShift in Cleng" << endl;
    bool success = true;

    int MX=Lat[0]->MX;
    int MY=Lat[0]->MY;
    int MZ=Lat[0]->MZ;


    int sum_of_elems = -5;
    int rand_part_index;
    vector<int> shift_XYZ = {0,0,0};

    rand_part_index = GetIntRandomValue(0, X.size()+1);
    cout << "rand_part_index:" << rand_part_index << endl;

//    while (sum_of_elems % 2 != 0 && any_of(shift_XYZ.begin(), shift_XYZ.end(), [&](int i){return i != 0;})) {
// TODO: Need here information about chainLenght and pathLenght ==> for generation new trials
// TODO: Need to add information about len of connections..

    while (sum_of_elems % 2 != 0 && count(shift_XYZ.begin(), shift_XYZ.end(), 0) != 1) {



        sum_of_elems = 0;

        for (int i = 0; i < 3; i++) {
            shift_XYZ[i] = GetIntRandomValue(-1, 1);
        }

        for_each(shift_XYZ.begin(), shift_XYZ.end(), [&] (int n) {
            sum_of_elems += n;
        });
    }
    cout << "Shift_inside:" << shift_XYZ[0] << " " << shift_XYZ[1] << " " << shift_XYZ[2] << endl;


    int changed = 0;
    for (int i=0; i < n_boxes; i++)
        if (changed != 1) {
            if (rand_part_index == P[2 * i]) {
                cout << "First" << endl;
                cout << " partx:" << X[P[2 * i]] << " " << X[P[2 * i]] + Sx[2 * i] * MX << endl;
                cout << " party:" << Y[P[2 * i]] << " " << Y[P[2 * i]] + Sy[2 * i] * MY << endl;
                cout << " partz:" << Z[P[2 * i]] << " " << Z[P[2 * i]] + Sz[2 * i] * MZ << endl;

                X[P[2 * i]] = X[P[2 * i]] + shift_XYZ[0];
                Y[P[2 * i]] = Y[P[2 * i]] + shift_XYZ[1];
                Z[P[2 * i]] = Z[P[2 * i]] + shift_XYZ[2];

                cout << "First_makeshift" << endl;
                cout << " partx:" << X[P[2 * i]] << " " << X[P[2 * i]] + Sx[2 * i] * MX << endl;
                cout << " party:" << Y[P[2 * i]] << " " << Y[P[2 * i]] + Sy[2 * i] * MY << endl;
                cout << " partz:" << Z[P[2 * i]] << " " << Z[P[2 * i]] + Sz[2 * i] * MZ << endl;

                changed += 1;
            }
            if (rand_part_index == P[2 * i + 1]) {
                cout << "Another" << endl;
                cout << " partx:" << X[P[2 * i + 1]] << " " << X[P[2 * i + 1]] + Sx[2 * i + 1] * MX << endl;
                cout << " party:" << Y[P[2 * i + 1]] << " " << Y[P[2 * i + 1]] + Sy[2 * i + 1] * MY << endl;
                cout << " partz:" << Z[P[2 * i + 1]] << " " << Z[P[2 * i + 1]] + Sz[2 * i + 1] * MZ << endl;

                X[P[2 * i + 1]] = X[P[2 * i + 1]] + shift_XYZ[0];
                Y[P[2 * i + 1]] = Y[P[2 * i + 1]] + shift_XYZ[1];
                Z[P[2 * i + 1]] = Z[P[2 * i + 1]] + shift_XYZ[2];

                cout << "Another makeshift" << endl;
                cout << " partx:" << X[P[2 * i + 1]] << " " << X[P[2 * i + 1]] + Sx[2 * i + 1] * MX << endl;
                cout << " party:" << Y[P[2 * i + 1]] << " " << Y[P[2 * i + 1]] + Sy[2 * i + 1] * MY << endl;
                cout << " partz:" << Z[P[2 * i + 1]] << " " << Z[P[2 * i + 1]] + Sz[2 * i + 1] * MZ << endl;

                changed += 1;
            }
        }

    cout << "Shift:" << shift_XYZ[0] << " " << shift_XYZ[1] << " " << shift_XYZ[2] << endl;

    return success;
}


bool Cleng::MonteCarlo() {
  if (debug) cout << "Monte Carlo in Cleng" << endl;
	bool success;

	Real free_energy;
	Real free_energy_t;
	Real my_rand = GetRealRandomValue(0, 1);

// init system outlook
    success = CP(to_cleng);   // to cleng
    New[0]->Solve(true);    // solving
    free_energy = Sys[0]-> FreeEnergy; // storing free_energy
    WriteOutput(0);  // writing output

    for (int i = 1; i < MCS; i++) { // loop for trials
        success=CP(to_cleng);
        MakeShift();
        success=CP(to_segment);
        New[0]->Solve(true);
        free_energy_t = Sys[0]-> FreeEnergy;
        if ((free_energy * my_rand) <= free_energy_t) {
            cout << "Accepted" << endl;
        } else {
            cout << "Deny" << endl;
        }
        WriteOutput(i);
    }
	return success;
}

void Cleng::PutParameter(string new_param) {
  KEYS.push_back(new_param);
}
string Cleng::GetValue(string parameter) {
  int i = 0;
  int length = PARAMETERS.size();
  while (i < length) {
    if (parameter == PARAMETERS[i]) {
      return VALUES[i];
    }
    i++;
  }
  return "";
}

void Cleng::PushOutput() {
	int* point;
	for (int i = 0; i < n_out; i++) {
		Out[i]->PointerVectorInt.clear();
		Out[i]->PointerVectorReal.clear();
		Out[i]->SizeVectorInt.clear();
		Out[i]->SizeVectorReal.clear();
  		Out[i]->strings.clear();
  		Out[i]->strings_value.clear();
 		Out[i]->bools.clear();
 	 	Out[i]->bools_value.clear();
 		Out[i]->Reals.clear();
 	 	Out[i]->Reals_value.clear();
  		Out[i]->ints.clear();
  		Out[i]->ints_value.clear();
		if (Out[i]->name=="ana" || Out[i]->name=="kal") Out[i]->push("t",t);
		if (Out[i]->name=="ana" || Out[i]->name=="kal") Out[i]->push("MCS",MCS);
		if (Out[i]->name=="ana" || Out[i]->name=="vec") { //example for putting an array of Reals of arbitrary length to output
			string s="vector;0"; //the keyword 'vector' is used for Reals; the value 0 is the first vector, use 1 for the next etc, 
			Out[i]->push("gn",s); //"gn" is the name that will appear in output file
			Out[i]->PointerVectorReal.push_back(Mol[clp_mol]->gn); //this is the pointer to the start of the 'vector' that is reported to output.
			Out[i]->SizeVectorReal.push_back(sizeof(Mol[clp_mol]->gn)); //this is the size of the 'vector' that is reported to output
		}
		if (Out[i]->name=="ana" || Out[i]->name=="pos") { //example for putting 3 array's of Integers of arbitrary length to output
			string s="array;0";
			Out[i]->push("X",s);
			point=X.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(X.size());
			s="array;1";
			Out[i]->push("Y",s);
			point = Y.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Y.size());
			s="array;2";
			Out[i]->push("Z",s);
			point=Z.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Z.size());
		}
	
	}
}
