#include "system.h"
#include "tools_host.h"
#include <algorithm>
#include <cmath>

System::System(const Input* In_, Lattice* Lat_, vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_, string name_)
{
	Seg = Seg_;
	Mol = Mol_;
	In = In_;
	name = name_;
	Sta = Sta_;
	Rea = Rea_;
	lat=Lat_;
	NAMICS_DBG( "Constructor for system " << endl);
	KEYS.push_back("initial_guess");
	KEYS.push_back("guess_inputfile");
	KEYS.push_back("write_initial_guess");
	KEYS.push_back("X");
	KEYS.push_back("E");

		charged=false;
	grad_epsilon = false;
	all_system=false;
	first_pass=true;
	neutralizer=-1;
}
System::~System()
{
	NAMICS_DBG( "Destructor for system " << endl);
	DeAllocateMemory();
}
void System:: DeAllocateMemory(void){
	NAMICS_DBG( "DeAllocateMemory in system " << endl);
	if (!all_system) return;
	free(H_GrandPotentialDensity);
	free(H_FreeEnergyDensity);
	free(H_alpha);
	if (charged)
	{
		free(H_q);
		free(H_psi);
	}
  free(phitot);
  free(TEMP);
  free(KSAM);
  free(CHI);

  if (charged) {
    free(EE);
    free(E);
    free(psiMask);
    free(eps);
  }

all_system=false;
}

void System::AllocateMemory()
{
	NAMICS_DBG( "AllocateMemory in system " << endl);
	DeAllocateMemory();
	int M = lat->M;
	H_GrandPotentialDensity = (Real *)malloc(M * sizeof(Real));
	std::fill(H_GrandPotentialDensity,H_GrandPotentialDensity+M,0);

	H_FreeEnergyDensity = (Real *)malloc(M * sizeof(Real));
	H_alpha = (Real *)malloc(M * sizeof(Real));
	if (charged)
	{
		H_q = (Real *)malloc(M * sizeof(Real));
		H_psi = (Real *)malloc(M * sizeof(Real));
	}

  phitot = (Real*)malloc(M * sizeof(Real));
  alpha = H_alpha;
  if (charged) {
    psi = H_psi;
    q = H_q;
    eps = (Real*)malloc(M * sizeof(Real));
    EE = (Real*)malloc(M * sizeof(Real));
     E = (Real*)malloc(M * sizeof(Real));
    psiMask = (Real*)malloc(M * sizeof(Real));
  }
  KSAM = (Real*)malloc(M * sizeof(Real));
  FreeEnergyDensity = H_FreeEnergyDensity;
  GrandPotentialDensity = H_GrandPotentialDensity;
  TEMP = (Real*)malloc(M * sizeof(Real));
  std::fill_n(KSAM, M, 0);
  if (charged) {
    std::fill_n(psi, M, 0);
    std::fill_n(EE, M, 0);
    std::fill_n(E, M, 0);
  }
	n_mol = In->MolList.size();
	lat->AllocateMemory();
	int n_mon = In->MonList.size();
	for (int i = 0; i < n_mon; ++i)
		Seg[i]->AllocateMemory();
	for (int i = 0; i < n_mol; ++i)
		Mol[i]->AllocateMemory();
	CheckChi_values(In->MonList.size()); //Here CHI matrix is allocated.
	all_system=true;
}

bool System::generate_mask()
{
	NAMICS_DBG( "generate_mask in system " << endl);
	int M = lat->M;
	bool success = true;
	FrozenList.clear();
	int length = In->MonList.size();
	for (int i = 0; i < length; i++)
	{

		if (Seg[i]->freedom == "frozen") {
			FrozenList.push_back(i);
		}
	}

	std::fill_n(KSAM, M, 0);

	length = FrozenList.size();
	for (int i = 0; i < length; ++i)
	{
		for (int __i = 0; __i < (M); ++__i) (KSAM)[__i] += (Seg[FrozenList[i]]->MASK)[__i];
	}

	for (int __i = 0; __i < (M); ++__i) (KSAM)[__i] = ((KSAM)[__i] == 0) ? 1 : 0;

	Real accessible_volume = 0;
	if (lat->gradients < 3)
	{
		for (int i = 0; i < M; i++)
		{
			accessible_volume += KSAM[i] * lat->L[i];
		}
	}
	else
	{
		for (int __i = 0; __i < (M); ++__i) (accessible_volume) += (KSAM)[__i];
	}

	lat->Accesible_volume=accessible_volume;

	return success;
}

bool System::PrepareForCalculations(bool first_time)
{
	NAMICS_DBG( "PrepareForCalculations in System " << endl);

	bool success = true;
	int M = lat->M;
		success = generate_mask();

	n_mol = In->MolList.size();
	success = lat->PrepareForCalculations();
	int n_mon = In->MonList.size();

	for (int i = 0; i < n_mon; i++)
	{
		success = Seg[i]->PrepareForCalculations(KSAM,first_time);
	}
	for (int i = 0; i < n_mol; i++)
	{
		success = Mol[i]->PrepareForCalculations(KSAM);
	}
	if (first_time) {
  		if (charged) {
    			int length = FrozenList.size();
    			std::fill_n(psiMask, M, 0);
    			fixedPsi0 = false;
    			for (int i = 0; i < length; ++i) {
      				if (Seg[FrozenList[i]]->fixedPsi0) {
        				fixedPsi0 = true;
        				for (int __i = 0; __i < (M); ++__i) (psiMask)[__i] += (Seg[FrozenList[i]]->MASK)[__i];
      				}
    			}
    			Real eps=Seg[0]->epsilon;
    			for (int i=1; i<n_mon; i++) {
				if (Seg[i]->epsilon != eps) grad_epsilon = true;
    			}
  		}
	}

  return success;
}

bool System::MakeItsLists(void) {
	bool changed=false;
	int length = In->MonList.size();
	SysMonList.clear();
	int ItMonListLength=ItMonList.size();
	ItMonList.clear();
	int ItStateListLength=ItStateList.size();
	ItStateList.clear();

	length = In->MolList.size();
	int statelength = In->StateList.size();
	int i = 0;
	while (i < length)
	{
		int j = 0;
		int LENGTH = Mol[i]->MolMonList.size();
			while (j < LENGTH)
			{
				if (!In->InSet(SysMonList, Mol[i]->MolMonList[j]))
				{
					SysMonList.push_back(Mol[i]->MolMonList[j]);
					if (Seg[Mol[i]->MolMonList[j]]->state_name.size() < 1 && IsUnique(Mol[i]->MolMonList[j], -1))
					{
						ItMonList.push_back(Mol[i]->MolMonList[j]);
					}
				}
				j++;
			}
		i++;
	}
	for (int j = 0; j < statelength; j++)
	{
		if (IsUnique(-1, j))
		{
			ItStateList.push_back(j);
		}
	}

	if ((ItMonListLength-ItMonList.size()==0 && ItStateListLength-ItStateList.size()==0) || ItMonListLength+ItStateListLength==0) changed = false; else changed = true;

	return changed;
}

bool System::CheckInput(int start_)
{
	NAMICS_DBG( "CheckInput for system " << endl);
	start=start_;
	bool success = true;
	bool solvent_found = false;
	solvent = -1; //value -1 means no solvent defined.
	Real phibulktot = 0;
	success = In->CheckParameters("sys", name, start, KEYS, PARAMETERS);
	if (success)
	{
		success = CheckChi_values(In->MonList.size());

			MakeItsLists();

		int length = In->MolList.size();
		int i = 0;
		while (i < length)
		{
			if (Mol[i]->freedom == "free")
				phibulktot += Mol[i]->phibulk;
			if (Mol[i]->freedom == "solvent")
			{
				solvent_found = true;
				solvent = i;
			}
			i++;
		}


		if (!solvent_found && In->MolList.size()==1) {
			if (Mol[0]->IsPinned()) {
				phibulktot=1; cout <<"WARNING: no solvent found. Expecting solvent free 'brush'" << endl;
			} else {
				cout <<"Error: No solvent molecule found. One of your molecules must have 'freedom' : 'solvent'! " << endl;
				success=false;
			}
		}
		else
		if (!solvent_found || (phibulktot > 0.99999999 && phibulktot < 1.0000000001))
		{
			cout << "In system '" + name + "' the 'solvent' was not found, while the volume fractions of the bulk do not add up to unity. " << endl;
			success = false;
		}
			bool has_charged_molecule = false;
			int molecule_count = In->MolList.size();
			for (int i = 0; i < molecule_count; i++) {
				if (Mol[i]->IsCharged()) {
					has_charged_molecule = true;
					break;
				}
		}
		if (has_charged_molecule)
		{
			charged = true;
			neutralizer = -1;
			bool neutralizer_needed = false;

			int length = In->MolList.size();
			for (int i = 0; i < length; i++)
				if (Mol[i]->freedom == "neutralizer")
					neutralizer = i;

			if (neutralizer < 0)
			{
				Real phibulk = 0;
				for (int i = 0; i < length; i++)
				{
					if (Mol[i]->freedom == "free")
						phibulk += Mol[i]->Charge() * Mol[i]->phibulk;
					if (Mol[i]->freedom == "restricted" && Mol[i]->IsCharged() && !Mol[i]->IsPinned())
						neutralizer_needed = true;
				}

				if (neutralizer_needed || abs(phibulk) > 1e-4)
				{
					cout << "Neutralizer needed because some un-pinned molecules have freedom 'restricted' and/or overall charge density of 'free' molecules in the bulk is not equal to zero " << endl;
					success = false;
				}
			}
		}

			vector<string> options;
			initial_guess = "previous_result";
		if (GetValue("initial_guess").size() > 0)
		{
			options.clear();
			options.push_back("previous_result");
			options.push_back("file");
			options.push_back("none");
			ParseString(GetValue("initial_guess"), initial_guess, options, " Info about 'initial_guess' rejected;");
			if (initial_guess == "file")
			{
				if (GetValue("guess_inputfile").size() > 0)
				{
					guess_inputfile = In->ResolvePath(GetValue("guess_inputfile"));
				}
				else
				{
					success = false;
					cout << " When 'initial_guess' is set to 'file', you need to supply 'guess_inputfile', but this entry is missing. Problem terminated " << endl;
				}
			}
		}
		write_initial_guess = false;
		if (GetValue("write_initial_guess").size() > 0) {
			ParseBool(GetValue("write_initial_guess"),
			          write_initial_guess,
			          " Info about 'write_initial_guess' rejected; default: 'false' used.");
		}
	}

	if (In->StateList.size() > 1)
	{
		int num_of_Seg_with_states = 0;
		int num_of_Eqns = In->ReactionList.size();
		int num_of_alphabulk_fixed = 0;
		int num_of_states = In->StateList.size();
		for (int k = 0; k < num_of_states; k++)
			if (Sta[k]->fixed)
				num_of_alphabulk_fixed++;
		int length = In->MonList.size();
		for (int k = 0; k < length; k++)
			if (Seg[k]->state_name.size() > 1)
				num_of_Seg_with_states++;
		if (num_of_Seg_with_states + num_of_Eqns + num_of_alphabulk_fixed != num_of_states)
		{
			cout << " num_of_Seg_with_states+num_of_Eqns+num_of_alphabulk_fixed !=num_of_states" << endl;
			if (num_of_alphabulk_fixed == 0)
			{
				cout << " Consider to define for one of the states an alphabulk value " << endl;
			}
			else
			{
				if (num_of_alphabulk_fixed > 1)
				{
					cout << " possibly you have specified too many alphabulk values for multiple states " << endl;
				}
				else
				{
					if ((num_of_Seg_with_states + num_of_Eqns + num_of_alphabulk_fixed > num_of_states))
					{
						cout << " Possibly you have defined too many equations ... " << endl;
					}
					else
					{
						cout << " Possibly you have defined too few equations ... " << endl;
					}
				}
			}
			success = false;
		}
	}

	if (GetValue("E").size()>0)
	{
		if ( !(GetValue("E")=="chi" || GetValue("E")=="Chi" || GetValue("E")=="CHI") ) {
			cout <<" Only the FH chi-interactions are implemented. Use 'sys : sysname : E : chi'" << endl;
			cout <<" Only chi-contributions to E are generated." << endl;
		}
	}
	if (GetValue("X").size() > 0)
	{
		XmolList.clear();
		XstateList_1.clear();
		XstateList_2.clear();
		Xn_1.clear();
		string s = GetValue("X");
		vector<string> sub;
		vector<string> SUB;
		In->split(s, '-', sub);
		if (sub[0] != "F")
		{
			cout << "X is the characteristic function specified by user." << endl;
			cout << "Example of how a characteristic function is defined:" << endl;
			cout << "Case 1: no internal states: 'F - molname_1  -molname_2 - ...' " << endl;
			cout << "Here molname_1 etc are names of molecules in the system. " << endl;
			cout << "Case 2: internal states : F - molname_1 -(statename_1,statename_2,#number) - ... " << endl;
			cout << "Note that in ... you can add as many molname's and (...,...,...) combinations as you wish" << endl;
			cout << "Here F = Helmholtz energy. " << endl;
			cout << "and the statement '-molname_i' implies that 'n_i times mu_i' is subtracted from F. " << endl;
			cout << "The (statename_i,statename_j,n) implies that 'n times theta_i times mu_j' is subtracted from F " << endl;
			cout << endl;
			cout << "Error found: The first item is not the expected 'F' " << endl;
			success = false;
		}
		else
		{
			int length_sub = sub.size();
			int length_mol = In->MolList.size();
			int length_state = In->StateList.size();
			for (int i = 1; i < length_sub; i++)
			{
				SUB.clear();
				In->split(sub[i], ',', SUB);
				if (SUB.size() == 1)
				{ //want to see mol name
					bool found = false;
					for (int k = 0; k < length_mol; k++)
					{
						if (SUB[0] == In->MolList[k])
						{
							XmolList.push_back(k);
							found = true;
						}
					}
					if (!found)
					{
						success = false;
						cout << "In characteristic function X, the entry '" + SUB[0] + "' is not a molecule name " << endl;
					}
				}
				else
				{ //want to see (statename1,statename2)
					if (SUB.size() != 3)
					{
						cout << "In characteristic function X, the entry '" + sub[i] + "' is not recognised as '(statename_1, statename_2,n_1 )' " << endl;
						success = false;
					}
					else
					{
						bool found_1 = false, found_2 = false;
						for (int k = 0; k < length_state; k++)
						{
							if (SUB[0].substr(1, SUB[0].length() - 1) == In->StateList[k])
							{
								found_1 = true;
								XstateList_1.push_back(k);
							}
							if (SUB[1] == In->StateList[k])
							{
								found_2 = true;
								XstateList_2.push_back(k);
							}
						}
						int sto = ParseInt(SUB[2].substr(0, SUB[2].length() - 1), -1);
						if (sto < 0)
						{
							success = false;
							cout << "In characteristic function X, the entry '" + sub[i] + "' does not include a positive integer at the third argument. " << endl;
						}
						else
							Xn_1.push_back(sto);
						if (!(found_1 && found_2))
						{
							cout << "In characteristic function X, the entry '" + sub[i] + "' is not coding for '(statename_1,statename_2,n_1)" << endl;
							if (!found_1)
								cout << "first state name " + SUB[0].substr(1, SUB[0].length() - 1) + " does not exist" << endl;
							if (!found_2)
								cout << "second state name " + SUB[1] + " does not exist" << endl;
							success = false;
						}
					}
				}
			}
		}
	}

	int length = In->MonList.size();

	std::array<int, 6> bc{};
	for (int i = 0; i < length; i++) {
		if (Seg[i]->freedom == "frozen") {
			if (Seg[i]->frozen_at_bound>-1) {
				if (Seg[i]->valence != 0) {
					cout <<"Currently it is not allowed to put a charged frozen segment in boundary. Put this frozen segment inside the system instead. " << endl;  success=false;
				}
				bc[Seg[i]->frozen_at_bound]++;
			}
		}
	}
	if (lat->BC[0]=="surface" && bc[0] ==0) {cout <<"Lonely 'surface'. Specify a segment with frozen_range including the lowerboundary in x" << endl; success=false;}
	if (lat->BC[1]=="surface" && bc[1] ==0) {cout <<"Lonely 'surface'. Specify a segment with frozen_range including the lowerboundary in y" << endl; success=false;}
	if (lat->BC[2]=="surface" && bc[2] ==0) {cout <<"Lonely 'surface'. Specify a segment with frozen_range including the lowerboundary in z" << endl; success=false;}
	if (lat->BC[3]=="surface" && bc[3] ==0) {cout <<"Lonely 'surface'. Specify a segment with frozen_range including the upperboundary in x" << endl; success=false;}
	if (lat->BC[4]=="surface" && bc[4] ==0) {cout <<"Lonely 'surface'. Specify a segment with frozen_range including the upperboundary in y" << endl; success=false;}
	if (lat->BC[5]=="surface" && bc[5] ==0) {cout <<"Lonely 'surface'. Specify a segment with frozen_range including the upperboundary in z" << endl; success=false;}
	if (lat->BC[0]=="surface" && bc[0] >1) {cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the lowerboundary in x" << endl; success=false;}
	if (lat->BC[1]=="surface" && bc[1] >1) {cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the lowerboundary in y" << endl; success=false;}
	if (lat->BC[2]=="surface" && bc[2] >1) {cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the lowerboundary in z" << endl; success=false;}
	if (lat->BC[3]=="surface" && bc[3] >1) {cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the upperboundary in x" << endl; success=false;}
	if (lat->BC[4]=="surface" && bc[4] >1) {cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the upperboundary in y" << endl; success=false;}
	if (lat->BC[5]=="surface" && bc[5] >1) {cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the upperboundary in z" << endl; success=false;}

	return success;
}

bool System::IsUnique(int Segnr_, int Statenr_)
{
	NAMICS_DBG( "System::IsUnique: Segnr = " << Segnr_ << " Statenr = " << Statenr_ << endl);
	bool is_unique = true;
	bool is_equal = true;
	int Segnr = Segnr_, Statenr = Statenr_;
	int length = 0;
	if (Segnr > -1)
		length = Seg[Segnr]->chi.size();
	else
		length = Sta[Statenr]->chi.size();
	int itmonlength = ItMonList.size();
	int itstatelength = ItStateList.size();
	if (Statenr < 0)
		for (int i = 0; i < itmonlength; i++)
		{
			if (is_unique)
			{
				for (int k = 0; k < length; k++)
					if (Seg[Segnr]->chi[k] != Seg[ItMonList[i]]->chi[k])
						is_equal = false;
				if (is_equal)
				{
					is_unique = false;
					Seg[Segnr]->unique = false;
					Seg[Segnr]->seg_nr_of_copy = ItMonList[i];
				}
				else
					is_equal = true;
			}
		}
	else
	{
		for (int i = 0; i < itmonlength; i++)
		{
			if (is_unique)
			{
				for (int k = 0; k < length; k++)
					if (Sta[Statenr]->chi[k] != Seg[ItMonList[i]]->chi[k])
						is_equal = false;
				if (is_equal)
				{
					is_unique = false;
					Sta[Statenr]->unique = false;
					Sta[Statenr]->seg_nr_of_copy = ItMonList[i];
				}
				else
					is_equal = true;
			}
		}
		for (int i = 0; i < itstatelength; i++)
		{
			if (is_unique)
			{
				for (int k = 0; k < length; k++)
					if (Sta[Statenr]->chi[k] != Sta[ItStateList[i]]->chi[k])
						is_equal = false;
				if (is_equal)
				{
					is_unique = false;
					Sta[Statenr]->unique = false;
					Sta[Statenr]->state_nr_of_copy = ItStateList[i];
				}
				else
					is_equal = true;
			}
		}
	}
	return is_unique;
}

void System::PutParameter(string new_param)
{
	NAMICS_DBG( "PutParameter for system " << endl);
	KEYS.push_back(new_param);
}

string System::GetValue(string parameter)
{
	NAMICS_DBG( "GetValue " + parameter + " for system " << endl);
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

void System::push(string s, Real X)
{
	NAMICS_DBG( "push (Real) for system " << endl);
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void System::push(string s, int X)
{
	NAMICS_DBG( "push (int) for system " << endl);
	ints.push_back(s);
	ints_value.push_back(X);
}
void System::push(string s, bool X)
{
	NAMICS_DBG( "push (bool) for system " << endl);
	bools.push_back(s);
	bools_value.push_back(X);
}
void System::push(string s, string X)
{
	NAMICS_DBG( "push (string) for system " << endl);
	strings.push_back(s);
	strings_value.push_back(X);
}
void System::PushOutput()
{
	NAMICS_DBG( "PushOutput for system " << endl);
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("e", e);
	push("k_B", k_B);
	push("eps0", eps0);
	push("temperature", T);
	push("free_energy", FreeEnergy);
	push("grand_potential", GrandPotential);
	push("start",start);
	if (lat->gradients==1) {
		push("Laplace_pressure",-GrandPotentialDensity[lat->fjc]);
	}
	if (lat->gradients==2) {
		if (lat->BC[4]=="surface") {
			push("Laplace_pressure",-GrandPotentialDensity[lat->P(2*lat->fjc,(lat->MY+lat->fjc)/2)]);
		} else {
			push("Laplace_pressure",-GrandPotentialDensity[lat->P(2*lat->fjc,lat->MY)]);
		}
	}
	if (GetValue("E").size() >0)
	{
		Real sumE=0;
		int length= In->MonList.size();
		for (int i=0; i<length; i++)
		for (int j=i+1; j<length; j++) {
			Real Eij=GetE(i,j);
			push("I_"+Seg[i]->name+"_"+Seg[j]->name,Eij);
			push("I_"+Seg[j]->name+"_"+Seg[i]->name,Eij);
			sumE+=Eij*Seg[i]->chi[j];
		}
		push("E",sumE);
	}
	int n_seg=In->MonList.size();
	for (int i=0; i<n_seg; i++)
	for (int j=0; j<n_seg; j++){
		push("chi_"+Seg[i]->name+"_"+Seg[j]->name,CHI[i * n_seg + j]);
	}
	Real X = 0;
	if (Xn_1.size() > 0 || XmolList.size() > 0)
	{
		X = FreeEnergy;
		int length_mol = XmolList.size();
		for (int i = 0; i < length_mol; i++)
		{
			X -= Mol[XmolList[i]]->n * Mol[XmolList[i]]->Mu;
		}
		int length_state = Xn_1.size();
		Real mu = -999;
		for (int i = 0; i < length_state; i++)
		{
			length_mol = In->MolList.size();
			for (int k = 0; k < length_mol; k++)
			{
				if (Mol[k]->chainlength == 1)
				{
					int seg = Mol[k]->MolMonList[0];
					if (Seg[seg]->ns > 1)
					{
						for (int j = 0; j < Seg[seg]->ns; j++)
						{
							if (Seg[seg]->state_name[j] == In->StateList[XstateList_2[i]])
							{
								mu = Mol[k]->mu_state[j];
							}
						}
					}
				}
			}
			if (mu == -999)
			{
				cout << "Failed to find chemical potential for state " + In->StateList[XstateList_2[i]] + ": (not a monomer?) In characteristic function X, mu is set to zero." << endl;
				mu = 0;
			}
			X -= Seg[Sta[XstateList_1[i]]->mon_nr]->state_theta[Sta[XstateList_1[i]]->state_nr] * Xn_1[i] * mu;
		}
		push("X", X);
		cout << " X  = " << X << endl;
	}
	if (solvent>-1) push("solvent", Mol[solvent]->name);
	string s = "profile;0";
	push("alpha", s);
	s = "profile;1";
	push("GrandPotentialDensity", s);
	push("grand_potential_density", s);
	s = "profile;2";
	push("FreeEnergyDensity", s);
	push("free_energy_density", s);
	s = "profile;6";
	push("phitot", s);

	if (charged)
	{
		push("Dpsi",psi[lat->M-1]-psi[0]);
		s = "profile;3";
		push("psi", s);
		s = "profile;4";
		push("q", s);
		s = "profile;5";
		push("eps", s);
	}
}

Real *System::GetPointer(string s, int &SIZE)
{
	NAMICS_DBG( "GetPointer for system " << endl);
	vector<string> sub;
	SIZE = lat->M;
	In->split(s, ';', sub);
	if (sub[1] == "0")
		return H_alpha;
	if (sub[1] == "1")
		return H_GrandPotentialDensity;
	if (sub[1] == "2")
		return H_FreeEnergyDensity;
	if (sub[1] == "3")
		return psi;
	if (sub[1] == "4")
		return q;
	if (sub[1] == "5")
		return eps;
	if (sub[1] == "6")
		return phitot;
	return NULL;
}
int *System::GetPointerInt(string s, int &SIZE)
{
	(void)SIZE;
	NAMICS_DBG( "GetPointerInt for system " << endl);
	vector<string> sub;
	In->split(s, ';', sub);
	if (sub[0] == "array")
	{ //set SIZE and return pointer of int array
	}
	return NULL;
}

int System::GetValue(string prop, int &int_result, Real &Real_result, string &string_result)
{
	NAMICS_DBG( "GetValue (long) for system " << endl);
	int length = ints.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == ints[i])
		{
			int_result = ints_value[i];
			return 1;
		}
	}
	length = Reals.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == Reals[i])
		{
			Real_result = Reals_value[i];
			return 2;
		}
	}
	length = bools.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == bools[i])
		{
			if (bools_value[i])
				string_result = "true";
			else
				string_result = "false";
			return 3;
		}
	}
	length = strings.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == strings[i])
		{
			string_result = strings_value[i];
			return 3;
		}
	}
	return 0;
}

bool System::CheckChi_values(int n_seg)
{
	NAMICS_DBG( "CheckChi_values for system " << endl);
	bool success = true;
	CHI = (Real *)malloc(n_seg * n_seg * sizeof(Real));
	for (int i = 0; i < n_seg; i++)
		for (int k = 0; k < n_seg; k++)
		{
			CHI[i * n_seg + k] = ParseReal(Seg[i]->GetValue("chi_" + Seg[k]->name), 123);
		}
	for (int i = 0; i < n_seg; i++)
		for (int k = 0; k < n_seg; k++)
			if (CHI[i * n_seg + k] == 123)
				CHI[i * n_seg + k] = CHI[k * n_seg + i];
	for (int i = 0; i < n_seg; i++)
		for (int k = 0; k < n_seg; k++)
			if (CHI[i * n_seg + k] == 123)
				CHI[i * n_seg + k] = 0;
	for (int i = 0; i < n_seg; i++)
	{
		if (CHI[i * n_seg + i] != 0)
		{
			cout << "CHI-values for 'same'-segments (e.g. CHI(x,x) ) should be zero. " << endl;
			success = false;
		}
	}
		for (int i = 0; i < n_seg; i++)
			for (int k = i + 1; k < n_seg; k++)
				if (CHI[i * n_seg + k] != CHI[k * n_seg + i])
				{
					cout << "CHI-value symmetry violated: chi(" << Seg[i]->name << "," << Seg[k]->name << ") is not equal to chi(" << Seg[k]->name << "," << Seg[i]->name << ")" << endl;
					success = false;
				}

		for (int i = 0; i < n_seg; i++) {
			string NAME=Seg[i]->GetValue("set_equal_to");
			if (NAME.size()>0) {
				int segnr = -1;
				for (int j=0; j<n_seg; j++) {
					if (Seg[j]->name ==NAME) segnr=j;
				}
				if (segnr<0 || segnr ==i) {
					if (segnr < 0) cout <<"In segment " << Seg[i]->name << " 'set_to_seg' is rejected because the segment " << NAME << " was not found" << endl;
					else cout <<"In segment " << Seg[i]->name << " 'set_to_seg' is rejected because the segment " << NAME << " can not copied from itself...." << endl;
				} else {
					Seg[i]->epsilon = Seg[segnr]->epsilon;
					for (int k=0; k<n_seg; k++) {CHI[i*n_seg+k]=CHI[segnr*n_seg+k]; CHI[k*n_seg+i]=CHI[i*n_seg+k]; }
					CHI[i*n_seg+segnr]=CHI[segnr*n_seg+i]=0;
				}
			}
		}

	int n_segments = In->MonList.size();
	int n_states = In->StateList.size();
	if (n_states == 1)
		n_states = 0;
	int n_chi = n_segments + n_states;


	for (int i = 0; i < n_segments; i++)
		for (int j = 0; j < n_segments; j++)
		{
			if (Seg[i]->chi[j] == -999 && Seg[j]->chi[i] == -999)
			{
				Seg[i]->chi[j] = Seg[j]->chi[i] = 0;
			}
			else
			{
				if (Seg[i]->chi[j] != -999 && Seg[j]->chi[i] == -999)
				{
					Seg[j]->chi[i] = Seg[i]->chi[j];
				}
				else
				{
					if (Seg[i]->chi[j] == -999 && Seg[j]->chi[i] != -999)
					{
						Seg[i]->chi[j] = Seg[j]->chi[i];
					}
					else
					{
						if (Seg[i]->chi[j] != Seg[j]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Seg[i]->name << "," << Seg[j]->name << ") != chi (" << Seg[j]->name << "," << Seg[i]->name << ")" << endl;
						}
					}
				}
			}
		}

	for (int i = n_segments; i < n_chi; i++)
		for (int j = n_segments; j < n_chi; j++)
		{
			if (Sta[i - n_segments]->chi[j] == -999 && Sta[j - n_segments]->chi[i] == -999)
			{
				Sta[i - n_segments]->chi[j] = Sta[j - n_segments]->chi[i] = Seg[Sta[i - n_segments]->mon_nr]->chi[Sta[j - n_segments]->mon_nr];
			}
			else
			{
				if (Sta[i - n_segments]->chi[j] != -999 && Sta[j - n_segments]->chi[i] == -999)
				{
					Sta[j - n_segments]->chi[i] = Sta[i - n_segments]->chi[j];
				}
				else
				{
					if (Sta[i - n_segments]->chi[j] == -999 && Sta[j - n_segments]->chi[i] != -999)
					{
						Sta[i - n_segments]->chi[j] = Sta[j - n_segments]->chi[i];
					}
					else
					{
						if (Sta[i - n_segments]->chi[j] != Sta[j - n_segments]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Sta[i - n_segments]->name << "," << Sta[j - n_segments]->name << ") != chi (" << Sta[j - n_segments]->name << "," << Sta[i - n_segments]->name << ")" << endl;
						}
					}
				}
			}
		}

	for (int i = 0; i < n_segments; i++)
		for (int j = n_segments; j < n_chi; j++)
		{
			if (Seg[i]->chi[j] == -999 && Sta[j - n_segments]->chi[i] == -999)
			{
				Seg[i]->chi[j] = Sta[j - n_segments]->chi[i] = Seg[i]->chi[Sta[j - n_segments]->mon_nr];
			}
			else
			{
				if (Seg[i]->chi[j] != -999 && Sta[j - n_segments]->chi[i] == -999)
				{
					Sta[j - n_segments]->chi[i] = Seg[i]->chi[j];
				}
				else
				{
					if (Seg[i]->chi[j] == -999 && Sta[j - n_segments]->chi[i] != -999)
					{
						Seg[i]->chi[j] = Sta[j - n_segments]->chi[i];
					}
					else
					{
						if (Seg[i]->chi[j] != Sta[j - n_segments]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Seg[i]->name << "," << Sta[j - n_segments]->name << ") != chi (" << Sta[j - n_segments]->name << "," << Seg[i]->name << ")" << endl;
						}
					}
				}
			}
		}

	for (int i = n_segments; i < n_chi; i++)
		for (int j = 0; j < n_segments; j++)
		{
			if (Sta[i - n_segments]->chi[j] == -999 && Seg[j]->chi[i] == -999)
			{
				Sta[i - n_segments]->chi[j] = Seg[j]->chi[i] = Seg[j]->chi[Sta[i - n_segments]->mon_nr];
			}
			else
			{
				if (Sta[i - n_segments]->chi[j] != -999 && Seg[j]->chi[i] == -999)
				{
					Seg[j]->chi[i] = Sta[i - n_segments]->chi[j];
				}
				else
				{
					if (Sta[i - n_segments]->chi[j] == -999 && Seg[j]->chi[i] != -999)
					{
						Sta[i - n_segments]->chi[j] = Seg[j]->chi[i];
					}
					else
					{
						if (Sta[i - n_segments]->chi[j] != Seg[j]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Sta[i - n_segments]->name << "," << Seg[j]->name << ") != chi (" << Seg[j]->name << "," << Sta[i - n_segments]->name << ")" << endl;
						}
					}
				}
			}
			if (Sta[i - n_segments]->mon_nr == j && Seg[j]->chi[i] != 0 && Seg[j]->chi[i] != -999)
			{
				success = false;
				cout << " chi between mon-type and one of its states is not allowed for chi(" << Sta[i - n_segments]->name << "," << Seg[j]->name << ")" << endl;
			}
		}

	return success;
}

void System::DoElectrostatics(Real *g, Real *x)
{
	int M = lat->M;
	int n_seg = In->MonList.size();
	std::fill_n(q, M, 0);
	std::fill_n(eps, M, 0);
	for (int i = 0; i < n_seg; i++)
	{
		if (Seg[i]->ns < 2)
		{
			for (int __i = 0; __i < (M); ++__i) (q)[__i] += (Seg[i]->valence) * (Seg[i]->phi)[__i];
		}
		lat->set_bounds(Seg[i]->phi);
		for (int __i = 0; __i < (M); ++__i) (eps)[__i] += (Seg[i]->epsilon) * (Seg[i]->phi)[__i];
	}
	int statelistlength = In->StateList.size();
	for (int i = 0; i < statelistlength; i++)
	{
		for (int __i = 0; __i < (M); ++__i) (q)[__i] += (Sta[i]->valence) * (Seg[Sta[i]->mon_nr]->phi_state + Sta[i]->state_nr * M)[__i];
	}

	for (int __i = 0; __i < (M); ++__i) (q)[__i] = ((phitot)[__i] != 0) ? ((q)[__i] / (phitot)[__i]) : 0;
	for (int __i = 0; __i < (M); ++__i) (eps)[__i] = ((phitot)[__i] != 0) ? ((eps)[__i] / (phitot)[__i]) : 0;
	std::copy_n(x, M, psi); std::copy_n(psi, M, g);
	lat->set_M_bounds(psi);
	if (fixedPsi0) {
		int length=FrozenList.size();
		for (int i=0; i<length; i++) {
			Seg[FrozenList[i]]->UpdateValence(g,psi,q,eps,grad_epsilon);
		}
	}
}

void System:: ComputePhis(Real* x,bool first_time, Real residual) {
	NAMICS_DBG("ComputPhis in  system " << endl);
	PutU(x);
	PrepareForCalculations(first_time);
	ComputePhis(residual);
}

bool System:: PutU(Real* xx) {
NAMICS_DBG("PutU in  Solve " << endl);
	int M=lat->M;
	int itmonlistlength=ItMonList.size();
	int itstatelistlength=ItStateList.size();
	int monlistlength =In->MonList.size();
	int statelistlength=In->StateList.size();
	int k=0;

	int itpos=(itmonlistlength+itstatelistlength)*M;
	Real valence;
	Real *u;
	bool success=true;
	const auto add_segment_contributions = [&](Real* field, Segment* seg, Real segment_valence) {
		for (int __i = 0; __i < (M); ++__i) (field)[__i] += (seg->u_ext)[__i];
		if (charged){
			for (int __i = 0; __i < (M); ++__i) (field)[__i] += (-1.0*seg->epsilon) * (EE)[__i];
			if (segment_valence !=0) {
				for (int __i = 0; __i < (M); ++__i) (field)[__i] += (segment_valence) * (psi)[__i];
			}
		}
	};

	if (charged) {
		std::copy_n(xx+itpos, M, psi);
		lat->UpdateEE(EE,psi,E);
	}


	for (int i=0; i<itmonlistlength; i++) {
		int IM=ItMonList[i];
		u=Seg[IM]->u;
		std::copy_n(xx+k*M, M, u);
		valence=Seg[IM]->valence;
		add_segment_contributions(u, Seg[IM], valence);
		for (int j=0; j<monlistlength; j++) {
			if (Seg[j]->seg_nr_of_copy==IM && Seg[j]->ns<2) {
				u=Seg[j]->u;
				std::copy_n(xx+k*M, M, u);
				valence=Seg[j]->valence;
				add_segment_contributions(u, Seg[j], valence);
			}

		}
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->seg_nr_of_copy==IM) {
				u=Seg[Sta[j]->mon_nr]->u+Sta[j]->state_nr*M;
				std::copy_n(xx+k*M, M, u);
				valence=Sta[j]->valence;
				add_segment_contributions(u, Seg[Sta[j]->mon_nr], valence);
			}
		}
		k++;
	}

	for (int i=0; i<itstatelistlength; i++) {
		int IS=ItStateList[i];
		u=Seg[Sta[IS]->mon_nr]->u+(Sta[IS]->state_nr)*M;
		std::copy_n(xx+k*M, M, u);
		valence=Sta[IS]->valence;
		add_segment_contributions(u, Seg[Sta[IS]->mon_nr], valence);
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->state_nr_of_copy==IS) {
				u=Seg[Sta[j]->mon_nr]->u+Sta[j]->state_nr*M;
				std::copy_n(xx+k*M, M, u);
				valence=Sta[j]->valence;
				add_segment_contributions(u, Seg[Sta[j]->mon_nr], valence);
			}
		}
		k++;
	}
	if (charged) itpos +=M;

	return success;
}

void System::Classical_residual(Real* x,Real*g,Real residual, int iterations, int iv){
NAMICS_DBG("Classical_residuals in scf mode in system " << endl);
	int M=lat->M;
	Real chi;
	int mon_length = In->MonList.size(); //also frozen segments
	int i,k;

	int itmonlistlength=ItMonList.size();
	int state_length = In->StateList.size();
	int itstatelistlength=ItStateList.size();

	std::copy_n(x, iv, g);
	ComputePhis(x,iterations==0,residual);
 	std::fill_n(alpha, M, 0);

	for (i=0; i<itmonlistlength; i++) {
		for (int __i = 0; __i < (M); ++__i) (g+i*M)[__i] += (Seg[ItMonList[i]]->u_ext)[__i];
		for (k=0; k<mon_length; k++) {
			if (Seg[k]->ns<2) {
				chi =Seg[ItMonList[i]]->chi[k];
				if (chi!=0) {
					for (int __i = 0; __i < (M); ++__i) if ((phitot)[__i] > 0) (g+i*M)[__i] = (g+i*M)[__i] - (chi) * (((Seg[k]->phi_side)[__i] / (phitot)[__i]) - (Seg[k]->phibulk));
				}
			}
		}
		for (k=0; k<state_length; k++) {
			chi =Seg[ItMonList[i]]->chi[mon_length+k];
			if (chi!=0) {
				for (int __i = 0; __i < (M); ++__i) if ((phitot)[__i] > 0) (g+i*M)[__i] = (g+i*M)[__i] - (chi) * (((Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr*M)[__i] / (phitot)[__i]) - (Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr]));
			}
		}
	}
	for (i=0; i<itmonlistlength; i++) for (int __i = 0; __i < (M); ++__i) (alpha)[__i] += (g+i*M)[__i];

	for (i=0; i<itstatelistlength; i++) {
		for (k=0; k<mon_length; k++) {
			if (Seg[k]->ns<2) {
				chi =Sta[ItStateList[i]]->chi[k];
				if (chi!=0) {
					for (int __i = 0; __i < (M); ++__i) if ((phitot)[__i] > 0) (g+(itmonlistlength+i)*M)[__i] = (g+(itmonlistlength+i)*M)[__i] - (chi) * (((Seg[k]->phi_side)[__i] / (phitot)[__i]) - (Seg[k]->phibulk));
				}
			}
		}


		for (k=0; k<state_length; k++) {
			chi =Sta[ItStateList[i]]->chi[mon_length+k];
			if (chi!=0) {
				for (int __i = 0; __i < (M); ++__i) if ((phitot)[__i] > 0) (g+(itmonlistlength+i)*M)[__i] = (g+(itmonlistlength+i)*M)[__i] - (chi) * (((Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr*M)[__i] / (phitot)[__i]) - (Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr]));
			}
		}
	}
	for (i=0; i<itstatelistlength; i++) for (int __i = 0; __i < (M); ++__i) (alpha)[__i] += (g+(itmonlistlength+i)*M)[__i];
	for (int __i = 0; __i < (M); ++__i) (alpha)[__i] *= (1.0/(itmonlistlength+itstatelistlength));
	for (i=0; i<itmonlistlength; i++) {
		for (int __i = 0; __i < (M); ++__i) (g+i*M)[__i] = ((phitot)[__i] > 0) ? ((g+i*M)[__i] - (alpha)[__i] + 1 / (phitot)[__i] - 1.0) : 0;
		lat->remove_bounds(g+i*M);
		for (int __i = 0; __i < (M); ++__i) (g+i*M)[__i] = (g+i*M)[__i] * (KSAM)[__i];
	}
	for (i=0; i<itstatelistlength; i++) {
		for (int __i = 0; __i < (M); ++__i) (g+(itmonlistlength+i)*M)[__i] = ((phitot)[__i] > 0) ? ((g+(itmonlistlength+i)*M)[__i] - (alpha)[__i] + 1 / (phitot)[__i] - 1.0) : 0;
		lat->remove_bounds(g+(itmonlistlength+i)*M);
		for (int __i = 0; __i < (M); ++__i) (g+(itmonlistlength+i)*M)[__i] = (g+(itmonlistlength+i)*M)[__i] * (KSAM)[__i];
	}

	int itpos=(itmonlistlength+itstatelistlength)*M;

	if (charged) {
		std::copy_n(x+itpos, M, g+itpos);
		DoElectrostatics(g+itpos,x+itpos);
		lat->set_M_bounds(psi);
		if (M > 1) psi[0] = psi[1];
		lat->UpdatePsi(g+itpos,psi,q,eps,psiMask,grad_epsilon,fixedPsi0);
		lat->remove_bounds(g+itpos);
		itpos+=M;
	}
}


bool System::ComputePhis(Real residual){
NAMICS_DBG("ComputePhis in system" << endl);
	int M= lat->M;
	Real A=0, B=0; //A should contain sum_phi*charge; B should contain sum_phi
	bool success=true;
	std::fill_n(phitot, M, 0);

	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		Real *phi_frozen=Seg[FrozenList[i]]->phi;
		for (int __i = 0; __i < (M); ++__i) (phitot)[__i] += (phi_frozen)[__i];
	}

	for (int i=0; i<n_mol; i++) {
		success = Mol[i]->ComputePhi();
	}

	for (int i = 0; i < n_mol; i++)
	{
		Real norm = 0;
		if (Mol[i]->freedom == "free")
		{
			norm = Mol[i]->phibulk / Mol[i]->chainlength;
			Mol[i]->n = norm * Mol[i]->GN;

			A += Mol[i]->phibulk * Mol[i]->Charge();
			B += Mol[i]->phibulk;
		}

		if (Mol[i]->freedom == "restricted")
		{
			if (Mol[i]->GN > 0)
			{
				norm = Mol[i]->n / Mol[i]->GN;
				if (Mol[i]->IsPinned())
				{
					Mol[i]->phibulk = 0;
				}
				else
				{
					Mol[i]->phibulk = Mol[i]->chainlength * norm;
					A += Mol[i]->phibulk * Mol[i]->Charge();
					B += Mol[i]->phibulk;
				}
			}
			else
				{
					norm = 0;
					cout << "GN for molecule " << i << " is not larger than zero..." << endl;
					cout << "Consider compiling with LongReal enabled if this is an overflow issue." << endl;
					throw - 1;
				}
		}

		if (Mol[i]->IsPinned())
		{
			if (Mol[i]->GN > 0){
				if (Mol[i]->n ==0) Mol[i]->n =1;
				norm = Mol[i]->n / Mol[i]->GN;
			} else
			{
				norm = 0;
				cout << "GN for molecule " << i << " is not larger than zero..." << endl;
			}
			Mol[i]->phibulk = 0;
		}
		int k = 0;
		Mol[i]->norm = norm;
		int length = Mol[i]->MolMonList.size();

		while (k < length)
		{
			if (Mol[i]->freedom != "frozen")
			{
				Real *phi = Mol[i]->phi + k * M;
				Real *G1 = Seg[Mol[i]->MolMonList[k]]->G1;
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] = ((G1)[__i] != 0) ? ((phi)[__i] / (G1)[__i]) : 0;
				if (norm > 0) {
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (norm);
				}

				if (debug)
				{
					Real sum;
					(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (phi)[__i];
					NAMICS_DBG("Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl);
				}
			}
			k++;
		}

	}
	if (charged && neutralizer > -1)
		{
		if (Mol[neutralizer]->Charge()==Mol[solvent]->Charge()) {
			cout << "WARNING: solvent charge equals neutralizer charge; outcome problematic...." << endl;
		} else Mol[neutralizer]->phibulk= ((B-1.0)*Mol[solvent]->Charge() -A)/(Mol[neutralizer]->Charge()-Mol[solvent]->Charge());

		if (Mol[neutralizer]->phibulk<0) {
			cout << "WARNING: neutralizer has negative phibulk. Consider changing neutralizer...: outcome problematic...." << endl;
cout <<"A is " << A << endl;
for (int j=0; j<n_mol; j++) {
	cout << " mol : " << Mol[j]->name << " phibulk " << Mol[j]->phibulk << endl;

}


		}
		B += Mol[neutralizer]->phibulk;
		Real norm = Mol[neutralizer]->phibulk / Mol[neutralizer]->chainlength;
		Mol[neutralizer]->n = norm * Mol[neutralizer]->GN;
		Mol[neutralizer]->theta = Mol[neutralizer]->n * Mol[neutralizer]->chainlength;
		Mol[neutralizer]->norm = norm;
	}

	if (solvent>-1) {
		Mol[solvent]->phibulk = 1.0 - B;
		if (Mol[solvent]->phibulk < 0)
		{
			cout << "WARNING: solvent has negative phibulk. outcome problematic " << endl;
			throw - 4;
		}

		Real norm = Mol[solvent]->phibulk / Mol[solvent]->chainlength;
		Mol[solvent]->n = norm * Mol[solvent]->GN;
		Mol[solvent]->theta = Mol[solvent]->n * Mol[solvent]->chainlength;
		Mol[solvent]->norm = norm;

		int k = 0;

		length = Mol[solvent]->MolMonList.size();
		while (k < length) {
			Real *phi = Mol[solvent]->phi + k * M;
			if (norm > 0) {
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (norm);
			}
			if (debug)
			{
				Real sum;
				(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (phi)[__i];
				NAMICS_DBG("Sumphi in mol " << solvent << "for mon " << k << ":" << sum << endl);
			}
			k++;
		}
	}
	if (charged && neutralizer > -1)
	{
		int k = 0;
		length = Mol[neutralizer]->MolMonList.size();
		while (k < length)
		{
			Real *phi = Mol[neutralizer]->phi + k * M;
			if (Mol[neutralizer]->norm > 0) {
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (Mol[neutralizer]->norm);
			}
			if (debug)
			{
				Real sum;
				(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (phi)[__i];
				NAMICS_DBG("Sumphi in mol " << neutralizer << "for mon " << k << ":" << sum << endl);
			}
			k++;
		}
	}

	for (int i = 0; i < n_mol; i++)
	{

		int length = Mol[i]->MolMonList.size();
		int k = 0;
		while (k < length)
		{
			Real *phi_mon = Seg[Mol[i]->MolMonList[k]]->phi;
			Real *mol_phitot = Mol[i]->phitot;
			Real *phi_molmon = Mol[i]->phi + k * M;
			for (int __i = 0; __i < (M); ++__i) (phi_mon)[__i] += (phi_molmon)[__i];
			for (int __i = 0; __i < (M); ++__i) (phitot)[__i] += (phi_molmon)[__i];
			for (int __i = 0; __i < (M); ++__i) (mol_phitot)[__i] += (phi_molmon)[__i];
			Seg[Mol[i]->MolMonList[k]]->phibulk += Mol[i]->fraction(Mol[i]->MolMonList[k]) * Mol[i]->phibulk;
			k++;
		}
		}

	int n_seg = In->MonList.size();
	for (int i = 0; i < n_seg; i++) {
		lat->set_bounds(Seg[i]->phi);
	}

	for (int i = 0; i < n_seg; i++) {
		Seg[i]->SetPhiSide();
	}
	return success;
}

bool System::CheckResults(bool e_info_)
{
	NAMICS_DBG( "CheckResults for system " << endl);

	bool e_info = e_info_;
	bool success = true;

	FreeEnergy = GetFreeEnergy();
	GrandPotential = GetGrandPotential();
	CreateMu(lat->M);
	int M = lat->M;
	for (int i = 0; i < n_mol; ++i) {
		if (!std::isfinite(Mol[i]->n)) {
			cerr << "Detected invalid n in computed solver state for molecule " << Mol[i]->name << "." << endl;
			return false;
		}
		if (!std::isfinite(Mol[i]->theta)) {
			cerr << "Detected invalid theta in computed solver state for molecule " << Mol[i]->name << "." << endl;
			return false;
		}
		if (!std::isfinite(Mol[i]->phibulk)) {
			cerr << "Detected invalid phibulk in computed solver state for molecule " << Mol[i]->name << "." << endl;
			return false;
		}
		if (!std::isfinite(Mol[i]->GN)) {
			cerr << "Detected invalid GN in computed solver state for molecule " << Mol[i]->name << "." << endl;
			return false;
		}
		for (int j = 0; j < M; ++j) {
			if (!std::isfinite(Mol[i]->phitot[j])) {
				cerr << "Detected invalid phitot in computed solver state for molecule " << Mol[i]->name << "." << endl;
				return false;
			}
		}
		const int n_molmon = Mol[i]->MolMonList.size();
		for (int j = 0; j < n_molmon * M; ++j) {
			if (!std::isfinite(Mol[i]->phi[j])) {
				cerr << "Detected invalid phi in computed solver state for molecule " << Mol[i]->name << "." << endl;
				return false;
			}
		}
	}

	if (e_info && first_pass)
	{
		cout << "free energy                 = " << FreeEnergy << endl;
		cout << "grand potential             = " << GrandPotential << endl;
	}
	Real n_times_mu = 0;
	for (int i = 0; i < n_mol; i++)
	{
		Real Mu = Mol[i]->Mu;
		Real n = Mol[i]->n;
		n_times_mu += n * Mu;
	}
	if (e_info && first_pass)
	{
		cout << "free energy     (GP + n*mu) = " << GrandPotential + n_times_mu << endl;
		cout << "grand potential (F - n*mu)  = " << FreeEnergy - n_times_mu << endl<<endl;;
		for (int i = 0; i < n_mol; i++)
		{
			int n_molmon = Mol[i]->MolMonList.size();
			Real theta_tot = Mol[i]->n * Mol[i]->chainlength;
			for (int j = 0; j < n_molmon; j++)
			{
				Real FRACTION = Mol[i]->fraction(Mol[i]->MolMonList[j]);
				Real THETA = lat->WeightedSum(Mol[i]->phi + j * M);
				cout << "MOL " << Mol[i]->name << " Fraction " << Seg[Mol[i]->MolMonList[j]]->name << ": " << FRACTION << "=?=" << THETA / theta_tot << " or " << THETA << " of " << theta_tot << endl;
			}
		}
	}
	first_pass=false;

	return success;
}

Real System::GetE(int Seg1, int Seg2)
{
	Real E=0;
	int M=lat->M;
        Real *temp = (Real *)malloc(M * sizeof(Real));
	Real *L=lat->L;
	Real *phi=Seg[Seg1]->phi;
	Real *side=Seg[Seg2]->phi_side;
	if (Seg1!=Seg2) {
		for (int __i = 0; __i < (M); ++__i) (temp)[__i] = (L)[__i] * (phi)[__i];
		for (int __i = 0; __i < (M); ++__i) (temp)[__i] = (temp)[__i] * (side)[__i];
		(E) = 0; for (int __i = 0; __i < (M); ++__i) (E) += (temp)[__i];
	}
	free (temp);
	return E;//only the contacts
}

Real System::GetFreeEnergy(void)
{ //eqn 2.91 of thesis of J.v.Male;
	NAMICS_DBG( "GetFreeEnergy for system " << endl);
	int M = lat->M;
	Real FreeEnergy = 0;
	Real *F = FreeEnergyDensity;
	Real constant = 0;
	int n_seg = In->MonList.size();
	int n_mol = In->MolList.size();
	int n_states = In->StateList.size();
	for (int i=0; i<n_mol; i++) {
		lat->remove_bounds(Mol[i]->phitot);
	}
	int n_mon = In->MonList.size();
	for (int i = 0; i < n_mon; i++)
	{
		if (Seg[i]->ns < 2)
			lat->remove_bounds(Seg[i]->phi_side);
		else
		{
			for (int j = 0; j < Seg[i]->ns; j++)
			{
				lat->remove_bounds(Seg[i]->phi_side + j * M);
			}
		}
	}

	std::fill_n(F, M, 0);

	for (int i = 0; i < n_mol; i++)
	{
		Real n = Mol[i]->n;
		Real GN = Mol[i]->GN;
		int N = Mol[i]->chainlength;
		Real *phi = Mol[i]->phitot;
		constant = log(N * n / GN) / N;
		std::copy_n(phi, M, TEMP);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (constant);
		for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
	}
	Real *phi;
	Real *phi_side;
	Real *g;
	Real chi;
	int n_sysmon = SysMonList.size();
	for (int j = 0; j < n_sysmon; j++)
	{
		int n_states=Seg[SysMonList[j]]->ns;

		if (n_states==1) {
			phi = Seg[SysMonList[j]]->phi;
			g = Seg[SysMonList[j]]->G1;
			for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = ((g)[__i] > 0) ? -log((g)[__i]) : 0;
			for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (TEMP)[__i];
			for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (-1);
			for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
		} else {
			for (int k=0; k<n_states; k++) {
				phi = Seg[SysMonList[j]]->phi_state + k* M;
				std::copy_n(Seg[SysMonList[j]]->u+k*M, M, TEMP);
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (TEMP)[__i];
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (-1);
				for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
			}
		}
	}

	for (int j = 0; j < n_seg; j++)
	{
		if (Seg[j]->ns < 2)
		{
			phi = Seg[j]->phi;
			for (int k = 0; k < n_seg; k++)
			{
				if (Seg[k]->ns < 2)
				{
					if (Seg[k]->freedom == "frozen")
						chi = Seg[j]->chi[k];
					else
						chi = Seg[j]->chi[k] / 2; //double counted.

					phi_side = Seg[k]->phi_side;
					if (!(Seg[j]->freedom == "frozen" || chi == 0))
					{
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
						for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
					}
				}
			}
			for (int i = 0; i < n_states; i++)
			{
					chi = Seg[j]->chi[n_seg + i]/2;

					phi_side = Seg[Sta[i]->mon_nr]->phi_side + Sta[i]->state_nr * M;
				if (!(Seg[j]->freedom == "frozen" || chi == 0))
				{
					for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
					for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
					for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
				}
			}
		}
	}



	for (int j = 0; j < n_states; j++)
	{
		phi = Seg[Sta[j]->mon_nr]->phi_state + Sta[j]->state_nr * M;
		for (int k = 0; k < n_seg; k++)
			if (Seg[k]->ns < 2)
			{
					chi = Sta[j]->chi[k];
					if (Seg[k]->freedom != "frozen") chi = chi / 2;


				phi_side = Seg[k]->phi_side;
				if (chi != 0)
				{
					for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
					for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
					for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
				}
			}
		for (int l = 0; l < n_states; l++)
		{
				chi = Sta[j]->chi[n_seg + l] / 2;

			phi_side = Seg[Sta[l]->mon_nr]->phi_side + Sta[l]->state_nr * M;
			if (chi != 0)
			{
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
				for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
			}
		}
	}

	for (int i = 0; i < n_mol; i++)
	{
		constant = 0;
		int n_molmon = Mol[i]->MolMonList.size();
		for (int j = 0; j < n_molmon; j++)
			for (int k = 0; k < n_molmon; k++)
			{
				Real fA = Mol[i]->fraction(Mol[i]->MolMonList[j]);
				Real fB = Mol[i]->fraction(Mol[i]->MolMonList[k]);
				Real chi = CHI[Mol[i]->MolMonList[j] * n_mon + Mol[i]->MolMonList[k]] / 2;
				constant -= fA * fB * chi;
			}
		Real *phi = Mol[i]->phitot;
		std::copy_n(phi, M, TEMP);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (constant);
		for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
	}
	for (int __i = 0; __i < (M); ++__i) (F)[__i] = (F)[__i] * (KSAM)[__i]; //clean up contributions in frozen sites.

std::fill_n(TEMP, M, 0);
	if (charged) {
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] += (q)[__i] * (psi)[__i];
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (0.5);
		for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
	}
	return FreeEnergy + lat->WeightedSum(F);
}

Real System::GetGrandPotential(void)
{ //Eqn 293
	NAMICS_DBG( "GetGrandPotential for system " << endl);
	int M = lat->M;
	Real *GP = GrandPotentialDensity;
	int n_mol = In->MolList.size();
	std::fill_n(GP, M, 0);

	for (int i = 0; i < n_mol; i++)
	{
		Real *phi = Mol[i]->phitot;
		Real phibulk = Mol[i]->phibulk;
		int N = Mol[i]->chainlength;
		std::copy_n(phi, M, TEMP);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (TEMP)[__i] + (-phibulk);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (1.0 / N); //GP has wrong sign. will be corrected at end of this routine;
		for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
	}

	for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (alpha)[__i];
	Real phibulkA;
	Real phibulkB;
	Real chi;
	Real *phi;
	Real *phi_side;
	Real *u_ext;
	int n_seg = In->MonList.size();
	int n_states = In->StateList.size();

	int n_mon = In->MonList.size();
	for (int i = 0; i < n_mon; i++) //if this is not done, in 3 gradients we have wrong results...
	{
		if (Seg[i]->ns < 2)
			lat->remove_bounds(Seg[i]->phi_side);
		else
		{
			for (int j = 0; j < Seg[i]->ns; j++)
			{
				lat->remove_bounds(Seg[i]->phi_side + j * M);
			}
		}
	}


		for (int j = 0; j < n_seg; j++)
			if (Seg[j]->freedom != "frozen")
			{
				if (Seg[j]->ns < 2)
				{
				phi = Seg[j]->phi;
				phibulkA = Seg[j]->phibulk;
					for (int k = 0; k < n_seg; k++)
					{
						if (Seg[k]->freedom != "frozen")
						{
						if (Seg[k]->ns < 2)
						{
								chi = Seg[j]->chi[k] / 2;

								phi_side = Seg[k]->phi_side;
							phibulkB = Seg[k]->phibulk;
							if (chi != 0)
							{
								for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
								for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (TEMP)[__i] + (-phibulkA * phibulkB);
								for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
								for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
							}
						}
					}
				}

				for (int k = 0; k < n_states; k++)
				{
						chi = Seg[j]->chi[n_seg + k] / 2;

						phi_side = Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr * M;
					phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
					if (chi != 0)
					{
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (TEMP)[__i] + (-phibulkA * phibulkB);
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
						for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
					}
				}
			}
		}
	for (int j=0; j< n_seg; j++)
	{
		phi=Seg[j]->phi;
		u_ext=Seg[j]->u_ext;
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (u_ext)[__i];
		for (int __i = 0; __i < (M); ++__i) (GP)[__i] -= (TEMP)[__i];
	}


	for (int j = 0; j < n_states; j++)
	{
		phi = Seg[Sta[j]->mon_nr]->phi_state + Sta[j]->state_nr * M;
		phibulkA = Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
			for (int k = 0; k < n_seg; k++)
				if (Seg[k]->freedom != "frozen")
			{
				if (Seg[k]->ns < 2)
				{
						chi = Sta[j]->chi[k] / 2;

						phibulkB = Seg[k]->phibulk;
					phi_side = Seg[k]->phi_side;
					if (chi != 0)
					{
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (TEMP)[__i] + (-phibulkA * phibulkB);
						for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
						for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
					}
				}
			}
		for (int k = 0; k < n_states; k++)
		{
				chi = Sta[j]->chi[n_seg + k] / 2;

			phi_side = Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr * M;
			phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
			if (chi != 0)
			{
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (TEMP)[__i] + (-phibulkA * phibulkB);
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
				for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
			}
		}
	}

	

	for (int __i = 0; __i < (M); ++__i) (GP)[__i] *= (-1.0); //correct the sign.

	std::fill_n(TEMP, M, 0);
if (charged) {
	for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (EE)[__i] * (eps)[__i];
	for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (-2.0);

	for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] += (q)[__i] * (psi)[__i];
	for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (-1.0/2.0);


		for (int __i = 0; __i < M; ++__i) GP[__i] += TEMP[__i];
		for (int __i = 0; __i < M; ++__i) GP[__i] = GP[__i] * KSAM[__i];
		for (int __i = 0; __i < M; ++__i) TEMP[__i] = q[__i] * KSAM[__i];
		for (int __i = 0; __i < M; ++__i) TEMP[__i] = q[__i] - TEMP[__i];
		for (int __i = 0; __i < M; ++__i) TEMP[__i] = TEMP[__i] * psi[__i];
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (0.5);
	for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
}

	if (!charged) for (int __i = 0; __i < (M); ++__i) (GP)[__i] = (GP)[__i] * (KSAM)[__i]; //exclude solid and frozen sites from GP.

	return  lat->WeightedSum(GP);

}

bool System::CreateMu(int pos)
{
	NAMICS_DBG( "CreateMu for system " << endl);
	int M=lat->M;
	bool success = true;
	Real constant;
	Real n;
	Real GN;
	int n_mol = In->MolList.size();
	int n_mon = In->MonList.size();
	for (int i = 0; i < n_mol; i++)
	{
		Real Mu = 0;
		Real NA = Mol[i]->chainlength;

			n = Mol[i]->n;
		GN = Mol[i]->GN;
		if (pos==M) Mu = log(NA * n / GN) + 1; else {
			Mu=log(Mol[i]->phitot[pos]) +1;
		}
		constant = 0;
		for (int k = 0; k < n_mol; k++)
		{
			Real NB = Mol[k]->chainlength;
				Real phibulkB;
		        if (pos==M) phibulkB=Mol[k]->phibulk; else phibulkB=Mol[k]->phitot[pos];
			constant += phibulkB / NB;
		}
		Mu = Mu - NA * constant;
		Real phibulkA;
		Real phibulkB;
		Real FA;
		Real FB;
		Real chi;
		int statelistlength = In->StateList.size();

		for (int j = 0; j < n_mon; j++)
		{
			if (Seg[j]->ns < 2)
			{
				if (pos==M) phibulkA = Seg[j]->phibulk; else phibulkA=Seg[j]->phi[pos];
					FA = Mol[i]->fraction(j);
				for (int k = 0; k < n_mon; k++)
				{
					if (Seg[k]->ns < 2)
					{
						if (pos==M) phibulkB = Seg[k]->phibulk; else phibulkB=Seg[k]->phi[pos];
							FB = Mol[i]->fraction(k);
						chi = Seg[j]->chi[k] / 2;
						Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
					}
				}
				for (int l = 0; l < statelistlength; l++)
				{     //include state contributions
					phibulkB = Seg[Sta[l]->mon_nr]->state_phibulk[Sta[l]->state_nr];
					FB = Mol[i]->fraction(Sta[l]->mon_nr) * Seg[Sta[l]->mon_nr]->state_alphabulk[Sta[l]->state_nr];
						chi = Seg[j]->chi[n_mon + l] / 2;
						Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);

				}
			}
		}
		for (int j = 0; j < statelistlength; j++)
		{ //adjust for steady state
			phibulkA = Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
			FA = Mol[i]->fraction(Sta[j]->mon_nr) * Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr];
				for (int k = 0; k < n_mon; k++)
					if (Seg[k]->ns < 2)
					{
						phibulkB = Seg[k]->phibulk;
						FB = Mol[i]->fraction(k);
						chi = Sta[j]->chi[k] / 2;
					Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
				}
			for (int k = 0; k < statelistlength; k++)
			{
				phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
					FB = Mol[i]->fraction(Sta[j]->mon_nr) * Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr];
					chi = Sta[j]->chi[n_mon + k] / 2;

				Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
			}
		}
 

		Mol[i]->Mu = Mu;
	}
	return success;
}
