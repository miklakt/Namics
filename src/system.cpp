#include "system.h"
#include "tools_host.h"
#include <algorithm>

System::System(const Input* In_, Lattice* Lat_, vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_, string name_)
{
	Seg = Seg_;
	Mol = Mol_;
	Lat = Lat_;
	In = In_;
	name = name_;
	Sta = Sta_;
	Rea = Rea_;
	lat=Lat;
	prepared = false;
	NAMICS_DBG( "Constructor for system " << endl);
  	KEYS.push_back("calculation_type");
	KEYS.push_back("constraint");
	KEYS.push_back("delta_range");
	KEYS.push_back("delta_range_units");
	KEYS.push_back("delta_inputfile");
	KEYS.push_back("delta_molecules");
	KEYS.push_back("phi_ratio");
	KEYS.push_back("initial_guess");
	KEYS.push_back("guess_inputfile");
	KEYS.push_back("final_guess");
	KEYS.push_back("guess_outputfile");
	KEYS.push_back("overflow_protection");
	KEYS.push_back("find_local_solution");
	KEYS.push_back("split");
	KEYS.push_back("X");
	KEYS.push_back("E");
	KEYS.push_back("compute_Gibbs_excess");
	KEYS.push_back("compute_kJ0");

	//  KEYS.push_back("guess-" + In->MonList[i]);
	charged=false;
	constraintfields=false;
  	boundaryless_volume=0;
	grad_epsilon = false;
	all_system=false;
	extra_constraints=0;
	local_solution=false;
	progress=0;
	old_residual = 10;
	do_blocks=false;
	first_pass=true;
	neutralizer=-1;
	pos_interface=0.0;
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
	if (constraintfields)
	{
		free(H_beta);
		free(H_BETA);
	}
  free(phitot);
  if (CalculationType=="steady_state") free(B_phitot);
  free(TEMP);
  free(KSAM);
  free(FILL);
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
	progress=0; old_residual=10;
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
	if (constraintfields)
	{
		H_beta = (Real *)malloc(M * sizeof(Real));
		std::fill(H_beta,H_beta+M,0);
		H_BETA = (Real *)malloc(M * sizeof(Real)); std::fill_n(H_BETA, M, 0);
		lat->FillMask(H_beta, px, py, pz, delta_inputfile);
	}

  phitot = (Real*)malloc(M * sizeof(Real));
  if (CalculationType=="steady_state") B_phitot = (Real*)malloc(M * sizeof(Real));
  alpha = H_alpha;
  if (charged) {
    psi = H_psi;
    q = H_q;
    eps = (Real*)malloc(M * sizeof(Real));
    EE = (Real*)malloc(M * sizeof(Real));
     E = (Real*)malloc(M * sizeof(Real));
    psiMask = (Real*)malloc(M * sizeof(Real));
  }
  if (constraintfields) {
	beta=H_beta;
	BETA=H_BETA;
  }
  KSAM = (Real*)malloc(M * sizeof(Real));
  FILL = (Real*)malloc(M * sizeof(Real));
  FreeEnergyDensity = H_FreeEnergyDensity;
  GrandPotentialDensity = H_GrandPotentialDensity;
  TEMP = (Real*)malloc(M * sizeof(Real));
  std::fill_n(KSAM, M, 0);
  std::fill_n(FILL, M, 0);
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
	extra_constraints=0;
	FrozenList.clear();
	int length = In->MonList.size();
	for (int i = 0; i < length; i++)
	{

		if (Seg[i]->freedom == "frozen") {
			FrozenList.push_back(i);
		}

		if (Seg[i]->constraints) extra_constraints+=Seg[i]->constraint_z.size();
	}

	// Tagged segments and frozen segments cannot occupy the same lattice sites.
	for (int ti = 0; ti < static_cast<int>(SysTagList.size()); ++ti) {
		const int tag_idx = SysTagList[ti];
		for (int fi = 0; fi < static_cast<int>(FrozenList.size()); ++fi) {
			const int frozen_idx = FrozenList[fi];
			for (int p = 0; p < M; ++p) {
				if (Seg[tag_idx]->MASK[p] > 0 && Seg[frozen_idx]->MASK[p] > 0) {
					cout << "Tagged segment '" << In->MonList[tag_idx]
					     << "' overlaps frozen segment '" << In->MonList[frozen_idx]
					     << "'. Please adjust tagged_range/frozen_range." << endl;
					return false;
				}
			}
		}
	}

	std::fill_n(KSAM, M, 0);

	length = FrozenList.size();
	for (int i = 0; i < length; ++i)
	{
		for (int __i = 0; __i < (M); ++__i) (KSAM)[__i] += (Seg[FrozenList[i]]->MASK)[__i];
	}

	length = SysTagList.size();
	for (int i = 0; i < length; ++i)
	{
		for (int __i = 0; __i < (M); ++__i) (KSAM)[__i] += (Seg[SysTagList[i]]->MASK)[__i];
	}

	length = SysClampList.size();
	for (int i = 0; i < length; ++i)
	{
		for (int __i = 0; __i < (M); ++__i) (KSAM)[__i] += (Seg[SysClampList[i]]->MASK)[__i];
	}


	for (int __i = 0; __i < (M); ++__i) (KSAM)[__i] = ((KSAM)[__i] == 0) ? 1 : 0;

	volume = 0;
	if (lat->gradients < 3)
	{
		for (int i = 0; i < M; i++)
		{
			volume += KSAM[i] * lat->L[i];
		}
	}
	else
	{
		for (int __i = 0; __i < (M); ++__i) (volume) += (KSAM)[__i];
	}

	// Boundary-free volume used by dynamic workflows.
	// This formula works across 1D/2D/3D lattices.
	this->boundaryless_volume = this->volume - ((2 * lat->gradients - 4) * lat->MX * lat->MY + 2 * lat->MX * lat->MZ + 2 * lat->MY * lat->MZ + (-2 + 2 * lat->gradients) * (lat->MX + lat->MY + lat->MZ) + pow(2, lat->gradients));

	lat->Accesible_volume=volume;

	return success;
}

bool System::PrepareForCalculations(bool first_time)
{
	NAMICS_DBG( "PrepareForCalculations in System " << endl);

	bool success = true;
	int M = lat->M;



		success = generate_mask();
		prepared = true;

	if (constraintfields)
	{
		for (int __i = 0; __i < (M); ++__i) (BETA)[__i] = exp(-(BETA)[__i]); // Beta wordt nu exp(-beta)
		//cin.get();
	}

	n_mol = In->MolList.size();
	success = lat->PrepareForCalculations();
	int n_mon = In->MonList.size();

	Filling=false;
	for (int i = 0; i < n_mol; i++){
		if (Mol[i]->Filling) {
			Filling=true;
			FillList.clear();
			std::fill_n(FILL, M, 0);
			Mol[i]->theta=0;
			Real frac=0;
			Real frac_0=0;
			Real pinned_v=0;
			int seg_pinned=Mol[i]->GetPinnedSeg();
			Real S1=0.0;
			Real S2=0.0;
			Mol[i]->FillRangesList.clear();
			for (int k=0; k<n_mon; k++) {
				S1=0.0;
				if (Seg[k]->freedom=="pinned") {
					for (int j=0; j<M; j++) S1+=Seg[k]->MASK[j]*Seg[seg_pinned]->MASK[j];
					S2=0.0; (S2) = 0; for (int __i = 0; __i < (M); ++__i) (S2) += (Seg[k]->MASK)[__i];
					if (S1==S2) Mol[i]->FillRangesList.push_back(k);
				}
			}
			int length=Mol[i]->FillRangesList.size();

			for (int j=0; j<length; j++) {
				frac=0;
				for (int __i = 0; __i < (M); ++__i) (FILL)[__i] += (Seg[Mol[i]->FillRangesList[j]]->MASK)[__i];
				FillList.push_back(Mol[i]->FillRangesList[j]);
				for (int k=0; k<n_mol;k++) {
					frac=Mol[k]->fraction(Mol[i]->FillRangesList[j]);
					if (k==i && frac>0) {
							if (frac_0==0) {
								frac_0=frac;
								pinned_v=Seg[Mol[i]->FillRangesList[j]]->PinnedVolume();
							} else {cout <<" Multiple pinned segments found in molecule that is filling the pinned_range" << endl; }
					} else Mol[i]->theta+=frac*Mol[k]->theta;
				}
			}
			if (frac_0==0) {
				success=false;
				cout <<"Error in computing theta for molecule " << Mol[i]->name << ". Possible pinned monomer of this molecule is not in the 'fill_range-of' list of monomers. " << endl;
			} else {
				Mol[i]->theta=(pinned_v-Mol[i]->theta)/frac_0;
				Mol[i]->n=Mol[i]->theta/Mol[i]->chainlength;
			success=false;
			}
			for (int i=0; i<M; i++) {
				if (FILL[i]>0) FILL[i]=1;
				lat->volume-=FILL[i];
			}

			for (int __i = 0; __i < (M); ++__i) (FILL)[__i] = ((FILL)[__i] == 0) ? 1 : 0;
		}
	}

	for (int i = 0; i < n_mon; i++)
	{
		success = Seg[i]->PrepareForCalculations(KSAM,first_time);
		if (Filling) {
			if (!In->InSet(FillList,i) && Seg[i]->freedom =="free") {
				for (int __i = 0; __i < (M); ++__i) (Seg[i]->G1)[__i] = (Seg[i]->G1)[__i] * (FILL)[__i];
			}
		}
	}
	for (int i = 0; i < n_mol; i++)
	{
		success = Mol[i]->PrepareForCalculations(KSAM);
	}
	if (first_time) {
		do_blocks=false;
		progress=0;
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
	StatelessMonList.clear();
	SysMolMonList.clear();
	SysMonList.clear();
	int ItMonListLength=ItMonList.size();
	ItMonList.clear();
	SysTagList.clear();
	SysClampList.clear();
	int ItStateListLength=ItStateList.size();
	ItStateList.clear();

	for (int i = 0; i < length; i++)
		if (Seg[i]->state_name.size() == 0)
			StatelessMonList.push_back(i);
	length = In->MolList.size();
	int statelength = In->StateList.size();
	int i = 0;
	while (i < length)
	{
		int j = 0;
		int LENGTH = Mol[i]->MolMonList.size();
		while (j < LENGTH)
		{
			SysMolMonList.push_back(Mol[i]->MolMonList[j]);
			if (!In->InSet(SysMonList, Mol[i]->MolMonList[j]))
			{
				if (Seg[Mol[i]->MolMonList[j]]->freedom != "tagged" && Seg[Mol[i]->MolMonList[j]]->freedom != "clamp")
				{
				SysMonList.push_back(Mol[i]->MolMonList[j]);
				if (Seg[Mol[i]->MolMonList[j]]->state_name.size() < 1 && IsUnique(Mol[i]->MolMonList[j], -1))
				{
					ItMonList.push_back(Mol[i]->MolMonList[j]);
					}
				}
			}
			if (Seg[Mol[i]->MolMonList[j]]->freedom == "tagged")
			{
				if (In->InSet(SysTagList, Mol[i]->MolMonList[j]))
				{
				}
				else
					SysTagList.push_back(Mol[i]->MolMonList[j]);
			}
			if (Seg[Mol[i]->MolMonList[j]]->freedom == "clamp")
				{
				if (In->InSet(SysClampList, Mol[i]->MolMonList[j]))
				{
				}
				else
					SysClampList.push_back(Mol[i]->MolMonList[j]);
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

string System::GetMonName(int mon_number_I)
{
	NAMICS_DBG( "GetMonName for system " << endl);
	return Seg[mon_number_I]->name;
}

bool System::CheckInput(int start_)
{
	NAMICS_DBG( "CheckInput for system " << endl);
	start=start_;
	bool success = true;
	bool solvent_found = false;
	tag_segment = -1;
	solvent = -1; //value -1 means no solvent defined. tag_segment=-1;
	Real phibulktot = 0;
	success = In->CheckParameters("sys", name, start, KEYS, PARAMETERS);
	if (success)
	{
		if (GetValue("find_local_solution").size()>0) {
			split = 2;
			local_solution=ParseBool(GetValue("find_local_solution"),false);
			if (local_solution) {
				if (lat->gradients!=3) {
					local_solution =false; cout << "find_local_solution is rejected as it requires 3 gradient system. " << endl;
					if (!(lat->MZ==2 || lat->MZ==4 || lat->MZ==8 || lat->MZ==16 ||lat->MZ==32 || lat->MZ==64 ||lat->MZ==128 || lat->MZ ==256)){
					      local_solution =false; cout<<"find_local_solution requires system size in z-direction equal to 2^n with n=1..8"<< endl;
					}
					if (!(lat->MY==2 || lat->MY==4 || lat->MY==8 || lat->MY==16 ||lat->MY==32 || lat->MY==64 ||lat->MY==128 || lat->MY ==256)){
					      local_solution =false; cout<<"find_local_solution requires system size in y-direction equal to 2^n with n=1..8"<< endl;
					}
					if (!(lat->MX==2 || lat->MX==4 || lat->MX==8 || lat->MX==16 ||lat->MX==32 || lat->MX==64 ||lat->MX==128 || lat->MX ==256)){
					      local_solution =false; cout<<"find_local_solution requires system size in x-direction equal to 2^n with n=1..8"<< endl;
					}
				}
				if (GetValue("split").size()>0) {
					split=ParseInt(GetValue("split"),2);
					if (!(split ==2 || split ==4 || split ==8 ||split ==16 || split==32 || split==64 || split ==128) ) {
						cout <<"Value for split should be 2^n, with n= 1,..,6. used split = 2 instead." << endl;
						split =2;
					}
					if (split > lat->MX || split > lat->MY || split > lat->MZ) {
						cout <<"Value for split can not exeed n_layers_x or n_layers_y or n_layers_z, value split=2 is used " << endl;
						split = 2;
					}

				}
			}
		}

		success = CheckChi_values(In->MonList.size());

#ifdef LongReal
		if (GetValue("overflow_protection").size()==0||GetValue("overflow_protection")=="false" || GetValue("overflow_protection")=="FALSE"){

			cout <<"The program is compiled for the use of 'long double' while 'overflow_protection' is not requested for;" << endl;
			cout <<"1. Turn on 'overflow_protection'." << endl;
			cout <<"2. Compile program without the #define 'LongReal' in namics.h. " << endl;
		}
#else
		if (GetValue("overflow_protection").size() > 0) {
			if (ParseBool(GetValue("overflow_protection"),true)) {
				cout<<"You request 'overflow_protection', but the program was not compiled with the #define LongReal" << endl;
				cout<<"1. Go to namics.h in the /src directory and turn on #define LongReal  ." <<endl;
				cout<<"2. Do not request 'overflow_protection'." << endl;
			}
		}

#endif

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
			if (Mol[i]->IsTagged())
			{
				tag_segment = Mol[i]->tag_segment;
				Mol[i]->n = 1.0 * Seg[tag_segment]->n_pos;
				Mol[i]->theta = Mol[i]->n * Mol[i]->chainlength;
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
		if (IsCharged())
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

		if (GetValue("constraint").size() > 0)
		{
			constraintfields = true;
			px.clear();
			py.clear();
			pz.clear();
			vector<string> constraints;
			constraints.push_back("delta");
			ConstraintType = "";
			if (!ParseString(GetValue("constraint"), ConstraintType, constraints, "Info about 'constraint' rejected"))
			{
				success = false;
			};
			if (ConstraintType == "delta")
			{
				//}
				if (GetValue("delta_range").size() > 0)
				{	int units=1;
					if (lat->fjc>1) {
						if (GetValue("delta_range_units").size() == 0) {
							cout <<"Because you have FJC-choices>3, you also need to specify the 'delta_range_units'. You can select 'bondlength' or 'gritsize'. " << endl;
							cout <<"Using bondlength units allows delta_range from 0 ... n_layers" << endl;
							cout <<"Using gritsize and e.g. FJC-choices=5 gives fjc=2 allows delta_range from fjc .. fjc (n_layers+1)-1, etc."  << endl; success=false;
						} else {
							vector<string> options;
							string bond_range_units;
							options.push_back("bondlength");
							options.push_back("gritsize");
							bond_range_units = GetValue("delta_range_units");
							if (bond_range_units=="bondlength") units = lat->fjc;
							else if (bond_range_units=="gritsize") units =1;
							else {
								cout << "Value for 'delta_range_units' not recognized. Use 'bondlength' or 'gritsize'. Depending on FJC-choices the delta_range can be larger for 'gritsize' than for 'bondlength'."<< endl;
								success=false; units =lat->fjc;
							}
						}
					} else {
 						units=1;
						if (GetValue("delta_range_units").size() > 0) {
							string delta_range_units=GetValue("delta_range_units");
							if (delta_range_units != "bondlength") cout << "Delta_range_units set to 'bondlength' because FJC-choices =3" << endl;
						}
					}
					string s = GetValue("delta_range");
					vector<string> sub;
					vector<string> set;
					vector<string> coor;
					In->split(s, ';', sub);
					int n_points = sub.size();
					for (int i = 0; i < n_points; i++)
					{
						set.clear();
						In->split(sub[i], '(', set);
						int length = set.size();
						if (length != 2)
						{
							if (length == 1 && set[0] == "file")
							{
								if (GetValue("delta_inputfile").size() > 0)
								{
									delta_inputfile = GetValue("delta_inputfile");
								}
								else
								{
									success = false;
									cout << "When 'delta_range' is set to 'file', you should provide a 'delta_inputfile'" << endl;
								}
							}
							else
							{
								success = false;
								cout << "In 'delta_range', for position " << i << "the expected '(x,y,z)' format was not found " << endl;
							}
						}
						else
						{
							coor.clear();
							In->split(set[1], ',', coor);
							int grad = lat->gradients;
							int corsize = coor.size();
							if (corsize != grad)
							{
								success = false;
								if (grad == 1)
									cout << "In 'delta_range', for position " << i << " the expected '(x)' format was not found " << endl;
								if (grad == 2)
									cout << "In 'delta_range', for position " << i << " the expected '(x,y)' format was not found " << endl;
								if (grad == 3)
									cout << "In 'delta_range', for position " << i << " the expected '(x,y,z)' format was not found " << endl;
							}
							else
							{
								int rr;
								rr=ParseInt(coor[0], -1)*units;
								if (rr<0 || rr>lat->MX) {cout << "Coordinate x for delta_range is out of bonds. " << endl; success=false; }
								else px.push_back(rr);
								if (grad > 1) {
									rr=ParseInt(coor[1], -1)*units;
									if (rr<0 || rr>lat->MY) {cout << "Coordinate y for delta_range is out of bonds. " << endl; success=false; }
									else py.push_back(rr);
								}
								if (grad > 2){
									rr=ParseInt(coor[2], -1)*units;
									if (rr<0 || rr>lat->MZ) {cout << "Coordinate z for delta_range is out of bonds. " << endl; success=false; }
									pz.push_back(rr);
								}
							}
						}
					}
				}
				else
				{
					success = false;
					cout << "When 'constraint' is set to 'delta', you should specify a 'delta_range' " << endl;
				}

				if (GetValue("delta_molecules").size() > 0)
				{
					string deltamols = GetValue("delta_molecules");
					vector<string> sub;
					In->split(deltamols, ';', sub);
					int length_sub = sub.size();
					if (length_sub != 2)
					{
						success = false;
						cout << " delta_molecules item should contain two 'molecule names' separated by a ';'" << endl;
					}
					else
					{
						DeltaMolList.clear();
						int length = In->MolList.size();
						for (int i = 0; i < length; i++)
						{
							if (sub[0] == Mol[i]->name)
								DeltaMolList.push_back(i);
							if (sub[1] == Mol[i]->name)
								DeltaMolList.push_back(i);
						}
						if (DeltaMolList.size() !=2) {success = false;
						cout << " In delta_molecules, two molecule names were expected but not found " << endl; return 0; }
						if (DeltaMolList[0] == DeltaMolList[1])
						{
							success = false;
							cout << " In delta_molecules you should specify two different names " << endl;
						}
						if (DeltaMolList.size() < 2)
						{
							success = false;
							cout << " In 'delta_molecules', one or more molecule names are not recognized" << endl;
						}
					}
				}
				else
				{
					success = false;
					cout << "When 'constraint' is set to 'delta', you should specify a set of 'delta_molecules' " << endl;
				}


				phi_ratio=-1.0;
				if(GetValue("phi_ratio").size()>0) {
					if (GetValue("phi_ratio")=="critical_ratio") {
						phi_ratio=1.0*Mol[DeltaMolList[0]]->chainlength/Mol[DeltaMolList[1]]->chainlength;
						if (phi_ratio>0) phi_ratio=sqrt(phi_ratio);
					}
					else phi_ratio=ParseReal(GetValue("phi_ratio"),-1);
					if (phi_ratio<0) {cout <<" phi_ratio shoud contain keyword 'critical_ratio' or a positive real number, typically 1. " << endl; success=false;}
				} else {
					success=false; cout <<"Please give a value for 'phi_ratio' (typically 1 or specify the keyword 'critical_ratio')" << endl;
				}
			}

		}


			vector<string> options;
			options.push_back("equilibrium");
			CalculationType = "equilibrium";
			if (GetValue("calculation_type").size() > 0)
			{
				if (!ParseString(GetValue("calculation_type"), CalculationType, options, " Info about calculation_type rejected; only 'equilibrium' is supported."))
				return false;
			}

			int num_of_gradient_settings=0;
			int num_of_mol = In->MolList.size();
			for (int i=0; i<num_of_mol; i++) {
				if (Mol[i]->freedom =="gradient") num_of_gradient_settings++;
			}

			if (num_of_gradient_settings>0) {
				cout << "Molecule freedom 'gradient' is not supported in this minimal build." << endl;
				return false;
			}

			if (CalculationType=="steady_state") {
				cout << "Calculation type 'steady_state' is not supported in this minimal build." << endl;
				return false;
			}

			if (false && CalculationType=="steady_state") {
				//Steady state is in development. For the time being this option is quite limited. In time some of these constraints will be lifted.
				if (lat->gradients>2) {
					cout <<"For 'calculation_type : steady_state' is currently limited to 1 gradient and 2 gradients calculations " << endl;
				return false;
			}
			if (lat->fjc != 1) {
				cout <<"For 'calculation_type : steady_state' the value for FJC_choices is limited to 3 " << endl;
				return false;
			}
			if (lat->BC[0] != "mirror") {
				cout <<"For 'calculation_type : steady_state' the setting for both 'lowerbound' and 'upperbound' must be 'mirror'. " << endl;
				return false;
			}
			for (int i=0; i<num_of_mol; i++) {
				int length=Mol[i]->MolMonList.size();
				for (int j=0; j<length; j++) {
					if (Seg[Mol[i]->MolMonList[j]]->used_in_mol_nr==-2) {
						Seg[Mol[i]->MolMonList[j]]->used_in_mol_nr=i;
						Seg[Mol[i]->MolMonList[j]]->B=Mol[i]->B;
						if (Mol[i]->phi_LB_X > 0) {
							Seg[Mol[i]->MolMonList[j]]->phi_LB_X=Mol[i]->phi_LB_X*Mol[i]->fraction(Mol[i]->MolMonList[j]);
							Seg[Mol[i]->MolMonList[j]]->phi_UB_X=Mol[i]->phi_UB_X*Mol[i]->fraction(Mol[i]->MolMonList[j]);
							Seg[Mol[i]->MolMonList[j]]->phi_LB_Y=Mol[i]->phi_LB_Y*Mol[i]->fraction(Mol[i]->MolMonList[j]);
							Seg[Mol[i]->MolMonList[j]]->phi_UB_Y=Mol[i]->phi_UB_Y*Mol[i]->fraction(Mol[i]->MolMonList[j]);
						}
					} else {
						cout <<"Multiple usage of seg " + Seg[Mol[i]->MolMonList[j]]->name + " in steady state system forbidden: each monomer may be used in only one molecule; it can not be used in multiple molecules. " << endl;
						return false;
					}
				}
			}
		}


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
					guess_inputfile = GetValue("guess_inputfile");
				}
				else
				{
					success = false;
					cout << " When 'initial_guess' is set to 'file', you need to supply 'guess_inputfile', but this entry is missing. Problem terminated " << endl;
				}
			}
		}
		final_guess = "next_problem";
		if (GetValue("final_guess").size() > 0)
		{
			options.clear();
			options.push_back("next_problem");
			options.push_back("file");
			if (!ParseString(GetValue("final_guess"), final_guess, options, " Info about 'final_guess' rejected; default: 'next_problem' used."))
			{
				final_guess = "next_problem";
			}
			if (final_guess == "file")
			{
				if (GetValue("guess_outputfile").size() > 0)
				{
					guess_outputfile = GetValue("guess_outputfile");
				}
				else
				{
					guess_outputfile = "";
					cout << "Filename not found for 'output_guess'. Default with inputfilename and extention '.outiv' is used. " << endl;
				}
			}
		}
	}

	internal_states = false;

	if (In->StateList.size() > 1)
	{
		internal_states = true;
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

	int *bc =(int*) malloc(6*sizeof(int)); std::fill(bc,bc+6,0);
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
	free(bc);

	return success;
}

bool System::IsUnique(int Segnr_, int Statenr_)
{
	NAMICS_DBG( "System::IsUnique: Segnr = " << Segnr_ << " Statenr = " << Statenr_ << endl);
	if (CalculationType=="steady_state") return true;
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

bool System::UpdateVarInfo(int step_nr) {
NAMICS_DBG("System:: UpdateVarInfo" << endl);
	bool success=true;
	switch(Var_scan_value) {
		case 0:
			if (scale=="exponential") {
				phi_ratio=pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value)+(1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				phi_ratio=Var_start_value+step_nr*Var_step;
			};
			cout <<"scanning ... sys : " + name + " : phi_ratio : " << phi_ratio << endl;
			break;

		default:
			cout <<"program error in System::UpdateVarInfo " << endl;
			break;
	}
	return success;
}

bool System::ResetInitValue() {
NAMICS_DBG("System:: ResetInitValue" << endl);
	bool success=true;
	cout <<"reset: ";
	switch (Var_scan_value) {
		case 0:
			phi_ratio=Var_start_value;
			cout <<"sys : " + name + " : phi_ratio : " << phi_ratio << endl;
			break;

		default:
			break;
	}
	return success;
}

int System::PutVarScan(Real step, Real end_value, int steps, string scale_) {
NAMICS_DBG("System:: PutVarScan" << endl);
	num_of_steps = -1;
	scale=scale_;
	Var_end_value=end_value;
	if (scale=="exponential") {
		Var_steps=steps; Var_step=0;
		if (steps==0) {
			cout <<"In var scan: the value of 'steps' is zero, this is not allowed" << endl; return -1;
		}
		if (Var_end_value*Var_start_value<0) {
			cout <<"In var scan: the product end_value*start_value <0. This is not allowed. " << endl; return -1;
		}
		if (Var_end_value > Var_start_value)
			num_of_steps=steps*log10(Var_end_value/Var_start_value);
		else
			num_of_steps=steps*log10(Var_start_value/Var_end_value);
	} else {
		Var_steps=0; Var_step=step;
		if (step==0) {
			cout <<"In var san: of system variable, the value of step can not be zero" << endl; return -1;
		}
		num_of_steps=(Var_end_value-Var_start_value)/step;

		if (num_of_steps<0) {
			cout <<"In var scan : (end_value-start_value)/step is negative. This is not allowed. Try changing the sign of 'step'. " << endl;
			return -1;
		}

	}
	return num_of_steps;
}

bool System::PutVarInfo(string Var_type_, string Var_target_, Real Var_target_value_)
{
	NAMICS_DBG( "System::PutVarInfo " << endl);
	bool success = true;
	if (Var_type_ =="scan") {
		Var_scan_value = -1;
		if (Var_target_=="phi_ratio") {
			Var_scan_value=0; Var_start_value=phi_ratio;
		}
		if (Var_scan_value<0) {
			success = false;
			cout << "Var scan " + Var_target_ + " rejected in PutVarInfo in System " << endl; return success;
		}
		return success;
	}

	Var_target = -1;
	if (Var_type_ != "target")
		success = false;
	if (Var_target_ == "free_energy")
		Var_target = 0;
	if (Var_target_ == "grand_potential")
		Var_target = 1;
	if (Var_target_ == "Laplace_pressure"){
			Var_target = 2;
	}
	if (Var_target < 0 || Var_target > 2 )
	{
		success = false;
		cout << "Var target " + Var_target_ + " rejected in PutVarInfo in System " << endl;
	}
	Var_target_value = Var_target_value_;
	if (Var_target_value < -1e4 || Var_target_value > 1e4) {
		success =false; cout <<"Var_target_value out of range " << endl;
	}
	return success;
}

Real System::GetError()
{
	NAMICS_DBG( "System::GetError " << endl);
	Real Error = 0;
	switch (Var_target)
	{
	case 0:
		Error = FreeEnergy - Var_target_value;
		break;
	case 1:
		Error = -1.0 * (GrandPotential - Var_target_value);
		break;
	case 2:
		Error = GrandPotentialDensity[lat->fjc]+GrandPotentialDensity[lat->M-2*lat->fjc]-Var_target_value;
		//cpush << " Error " << Error << endl;
		break;
	default:
		cout << "Program error in GetVarError" << endl;
		break;
	}
	return Error;
}

bool System::IsCharged()
{
	NAMICS_DBG( "System::IsCharged " << endl);
	bool success = false;
	int length = In->MolList.size();
	for (int i = 0; i < length; i++)
	{
		if (Mol[i]->IsCharged())
			success = true;
	}
	return success;
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


	if (GetValue("delta_range").size()>0) push("delta_range",GetValue("delta_range"));
	if (GetValue("phi_ratio").size()>0) push("phi_ratio",phi_ratio);
	int n_seg=In->MonList.size();
	for (int i=0; i<n_seg; i++)
	for (int j=0; j<n_seg; j++){
		push("chi_"+Seg[i]->name+"_"+Seg[j]->name,CHI[i * n_seg + j]);
	}
	if (GetValue("compute_kJ0").size()>0){
		int M=lat->M;
		if (lat->gradients==1 && lat->geometry=="planar") {
			if (pos_interface==0) pos_interface=M/2+0.5;
			push("kJ0", -lat->MomentPlanar(GrandPotentialDensity,1,pos_interface)/lat->fjc);
			pos_interface=0;
		} else {
			cout <<" 'compute_kJ0' requested but 'compute_kJ0' rejected because either geomety is not planar, or gradients = 1 or 'delta_range' not found " << endl;
		}

		if (lat->gradients == 1 && lat->geometry == "planar") {
			push("kbar", lat->MomentPlanar(GrandPotentialDensity,2,M/2+0.5)/pow(lat->fjc,2));
		}
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
	push("calculation_type", CalculationType);

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
	int n_mol = In->MolList.size();
	Real Sprod=0;
	for (int i=0; i<n_mol; i++) {
		Sprod += Mol[i]->J*Mol[i]->Delta_MU;
	}
	push("Sprod",Sprod);

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

int System::GetMonNr(string MonName) {
	int nr=-1;
	int length=In->MonList.size();
	for (int i=0; i<length; i++) {
		if (Seg[i]->name ==MonName) nr=i;
	}
	return nr;
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
		string NAME=Seg[i]->GetOriginal();
		if (NAME.size()>0) {
			int segnr = GetMonNr(NAME);
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

	//}



	int n_segments = In->MonList.size();
	int n_states = In->StateList.size();
	if (n_states == 1)
		n_states = 0;
	int n_chi = n_segments + n_states;

	//{for (int k=0; k<n_chi; k++) cout <<Seg[i]-> chi[k] << " "; cout << endl; }
	//{for (int k=0; k<n_chi; k++) cout <<Sta[i]-> chi[k] << " "; cout << endl; }

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

bool System:: Put_U(Real* xx){
	NAMICS_DBG( "Put_U in System" << endl);
	bool success=true;
	int M=lat->M;
	int itmonlistlength=ItMonList.size();
	for (int i=0; i<itmonlistlength; i++) {
		int IM=ItMonList[i];
	 	Real *u=Seg[IM]->u;
		std::copy_n(u, M, xx+i*M);
	}
	return success;
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

	if (charged) {
		std::copy_n(xx+itpos, M, psi);
		lat->UpdateEE(EE,psi,E);
	}


	for (int i=0; i<itmonlistlength; i++) {
		int IM=ItMonList[i];
		u=Seg[IM]->u;
		std::copy_n(xx+k*M, M, u);
		if (charged){
			for (int __i = 0; __i < (M); ++__i) (u)[__i] += (-1.0*Seg[IM]->epsilon) * (EE)[__i];
			valence=Seg[IM]->valence;
			if (valence !=0)
				for (int __i = 0; __i < (M); ++__i) (u)[__i] += (valence) * (psi)[__i];
		}
		for (int j=0; j<monlistlength; j++) {
			if (Seg[j]->seg_nr_of_copy==IM && Seg[j]->ns<2) {
				u=Seg[j]->u;
				std::copy_n(xx+k*M, M, u);
				if (charged){
					for (int __i = 0; __i < (M); ++__i) (u)[__i] += (-1.0*Seg[j]->epsilon) * (EE)[__i];
					valence=Seg[j]->valence;
					if (valence !=0)
						for (int __i = 0; __i < (M); ++__i) (u)[__i] += (valence) * (psi)[__i];
				}
			}

		}
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->seg_nr_of_copy==IM) {
				u=Seg[Sta[j]->mon_nr]->u+Sta[j]->state_nr*M;
				std::copy_n(xx+k*M, M, u);
				if (charged){
					for (int __i = 0; __i < (M); ++__i) (u)[__i] += (-1.0*Seg[Sta[j]->mon_nr]->epsilon) * (EE)[__i];
					valence=Sta[j]->valence;
					if (valence !=0)
						for (int __i = 0; __i < (M); ++__i) (u)[__i] += (valence) * (psi)[__i];


				}
			}
		}
		k++;
	}

	for (int i=0; i<itstatelistlength; i++) {
		int IS=ItStateList[i];
		u=Seg[Sta[IS]->mon_nr]->u+(Sta[IS]->state_nr)*M;
		std::copy_n(xx+k*M, M, u);
		if (charged){
			for (int __i = 0; __i < (M); ++__i) (u)[__i] += (-1.0*Seg[Sta[IS]->mon_nr]->epsilon) * (EE)[__i];
			valence=Sta[IS]->valence;
			if (valence !=0)
				for (int __i = 0; __i < (M); ++__i) (u)[__i] += (valence) * (psi)[__i];
		}
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->state_nr_of_copy==IS) {
				u=Seg[Sta[j]->mon_nr]->u+Sta[j]->state_nr*M;
				std::copy_n(xx+k*M, M, u);
				if (charged){
					for (int __i = 0; __i < (M); ++__i) (u)[__i] += (-1.0*Seg[Sta[j]->mon_nr]->epsilon) * (EE)[__i];
					valence=Sta[j]->valence;
					if (valence !=0)
						for (int __i = 0; __i < (M); ++__i) (u)[__i] += (valence) * (psi)[__i];
				}
			}
		}
		k++;
	}
	if (charged) itpos +=M;
	if (constraintfields) {std::copy_n(xx+itpos, M, BETA); itpos+=M;}
	if (extra_constraints>0) {
		int length = In->MonList.size();
		for (int i = 0; i < length; i++)
		{
			int constraint_size=Seg[i]->constraint_z.size();
			for (int k=0; k<constraint_size; k++) {
				itpos++;
				Seg[i]->Put_beta(k,xx[itpos-1]);
			}
		}
	}

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

	if (constraintfields) { //only works for two components...
		std::copy_n(Mol[DeltaMolList[1]]->phitot, M, g+itpos);
		for (int __i = 0; __i < (M); ++__i) (g+itpos)[__i] = (g+itpos)[__i] - (Mol[DeltaMolList[0]]->phitot)[__i];
		Real R = (phi_ratio-1)/(phi_ratio+1);
		for (int __i = 0; __i < (M); ++__i) (g+itpos)[__i] = (g+itpos)[__i] + (R);
		for (int __i = 0; __i < (M); ++__i) (g+itpos)[__i] = (g+itpos)[__i] * (beta)[__i];
		itpos+=M;
	}

	if (extra_constraints>0) {
		int length = In->MonList.size();
		for (int i = 0; i < length; i++)
		{
			int constraint_size=Seg[i]->constraint_z.size();
			for (int k=0; k<constraint_size; k++) {
				itpos++;
				g[itpos-1]=Seg[i]->Get_g(k);
			}
		}
	}
}


void System::Steady_residual(Real* x,Real*g,Real residual, int iterations, int iv){
NAMICS_DBG("steady_residuals in scf mode in system " << endl);
	int M=lat->M;
	Real chi;
	int mon_length = In->MonList.size(); //also frozen segments
	int i,k;
	//dphidt*=(Real*) malloc(iv*sizeof(Real);

	int itmonlistlength=ItMonList.size();
	int state_length = In->StateList.size();
	int itstatelistlength=ItStateList.size();

	if (itstatelistlength>0) cout <<"currently, internal states of segments incompatible with steady state " << endl;
	std::copy_n(x, iv, g);
	ComputePhis(x,iterations==0,residual);

	//	Seg[ItMonList[i]]->SetPhiSide();
	//}
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
	for (i=0; i<itmonlistlength; i++) std::copy_n(g+i*M, M, Seg[ItMonList[i]]->ALPHA);

	std::fill_n(g, iv, 0);
	Real Jtot=0;
	Segment* Seg0=Seg[ItMonList[0]];



	int gradients=lat->gradients;
	int MX=lat->MX;
	int MY=lat->MY;
	int JX=lat->JX;
	for (int z=1; z<M-1; z++) g[z]=1.0/phitot[z]-1.0;

	switch (gradients) {
		case 1:
			g[1]=Seg0->phi[0]/Seg0->phi[1]-1.0;
			//g[M-2]=Seg0->phi[M-1]/Seg0->phi[M-2]-1.0;
			g[MX]=Seg0->phi[MX+1]/Seg0->phi[MX]-1.0;
			break;
		case 2:
			for (int x=1; x<MX+1; x++) {
				g[x*JX+1]=Seg0->phi[x*JX+0]/Seg0->phi[x*JX+1]-1.0;
				g[x*JX+MY]=Seg0->phi[x*JX+MY+1]/Seg0->phi[x*JX+MY]-1.0;
			}
			for (int y=1; y<MY+1; y++) {
				g[JX+y]=Seg0->phi[0+y]/Seg0->phi[JX+y]-1.0;
				g[MX*JX+y]=Seg0->phi[(MX+1)*JX+y]/Seg0->phi[MX*JX+y]-1.0;
			}
			break;
		default:
			break;
	}
	for (int i =1 ; i<itmonlistlength; i++) {
		Segment* Segi=Seg[ItMonList[i]];
		Segi->J=0;
		for (int k =0; k<itmonlistlength; k++) {
			Segment* Segk=Seg[ItMonList[k]];
			if (i !=k) Segi->J += lat->DphiDt(g+i*M,B_phitot,Segi->phi,Segk->phi,Segi->ALPHA,Segk->ALPHA,Segi->B,Segk->B);
		}
		Jtot +=Segi->J;
	}
	Seg[ItMonList[0]]->J=-Jtot;
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
	if (constraintfields) {
		std::copy_n(Mol[DeltaMolList[1]]->phitot, M, g+itpos);
		for (int __i = 0; __i < (M); ++__i) (g+itpos)[__i] = (g+itpos)[__i] - (Mol[DeltaMolList[0]]->phitot)[__i];
		Real R = (phi_ratio-1)/(phi_ratio+1);
		for (int __i = 0; __i < (M); ++__i) (g+itpos)[__i] = (g+itpos)[__i] + (R);
		for (int __i = 0; __i < (M); ++__i) (g+itpos)[__i] = (g+itpos)[__i] * (beta)[__i];
		itpos+=M;
	}
	if (extra_constraints>0) {
		int length = In->MonList.size();
		for (int i = 0; i < length; i++)
		{
			int constraint_size=Seg[i]->constraint_z.size();
			for (int k=0; k<constraint_size; k++) {
				itpos++;
				g[itpos-1]=Seg[i]->Get_g(k);
			}
		}
	}
}

bool System::ComputePhis(Real residual){
NAMICS_DBG("ComputePhis in system" << endl);
	bool prepare_for_blocks=false;
	int M= lat->M;
	Real A=0, B=0; //A should contain sum_phi*charge; B should contain sum_phi
	bool success=true;
	std::fill_n(phitot, M, 0);

	if (local_solution) {
		if (residual/old_residual < 0.9) progress--; else progress++;
		if (progress < 0 ) progress=0;
		if (progress > 100 ) {
			progress=0;
			if (!do_blocks) {
			prepare_for_blocks=true;
			cout <<"no progress: trying fo find local solution" << endl;
			}
		}
		old_residual = residual;
	}

	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		Real *phi_frozen=Seg[FrozenList[i]]->phi;
		for (int __i = 0; __i < (M); ++__i) (phitot)[__i] += (phi_frozen)[__i];
	}

	for (int i=0; i<n_mol; i++) {
		if (constraintfields) {
			if (i==DeltaMolList[0]) {
				success=Mol[i]->ComputePhi(BETA,1);
			} else {
				if (i==DeltaMolList[1]) {
					success=Mol[i]->ComputePhi(BETA,-1);
				} else {
					success=Mol[i]->ComputePhi(BETA,0);
				}
			}
		}
		else {
			success = Mol[i]->ComputePhi(BETA, 0);
		}
	}

	for (int i = 0; i < n_mol; i++)
	{
		Real norm = 0;
		if (Mol[i]->freedom == "free" || Mol[i]->freedom == "gradient")
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
				cout << "Consider to turn on the overflow_protection (go to namics.h and #define LongReal, recompile using 'make'.)" << endl;
				throw - 1;
			}
		}

		if (Mol[i]->IsTagged() || Mol[i]->IsPinned())
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
		if (Mol[i]->IsClamped())
		{
			norm = 1;
			Mol[i]->phibulk = 0;
		}

		int k = 0;
		Mol[i]->norm = norm;
		int length = Mol[i]->MolMonList.size();

		while (k < length)
		{
			if (!(Seg[Mol[i]->MolMonList[k]]->freedom == "clamp" || Mol[i]->freedom == "frozen"))
			{
				Real *phi = Mol[i]->phi + k * M;
				Real *G1 = Seg[Mol[i]->MolMonList[k]]->G1;
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] = ((G1)[__i] != 0) ? ((phi)[__i] / (G1)[__i]) : 0;
				if (norm > 0)
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (norm);

					if (debug)
					{
						Real sum;
						(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (phi)[__i];
						NAMICS_DBG("Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl);
					}
			}
			k++;
		}

		if (Mol[i]->freedom == "range_restricted")
		{
			Real *phit = Mol[i]->phitot;
			std::fill_n(phit, M, 0);
			int k = 0;
			while (k < length)
			{
				Real *phi = Mol[i]->phi + k * M;
				for (int __i = 0; __i < (M); ++__i) (phit)[__i] += (phi)[__i];
				k++;
			}
			OverwriteA(phit, Mol[i]->R_mask, phit, M);
			//lat->remove_bounds(phit);
			Real theta = lat->ComputeTheta(phit);
			norm = Mol[i]->theta_range / theta;
			Mol[i]->norm = norm;
			Mol[i]->phibulk = Mol[i]->chainlength * norm;

			A += Mol[i]->phibulk * Mol[i]->Charge();
			B += Mol[i]->phibulk;

			Mol[i]->n = norm * Mol[i]->GN;
			Mol[i]->theta = Mol[i]->n * Mol[i]->chainlength;
			std::fill_n(phit, M, 0);
			k = 0;
			while (k < length)
			{
				Real *phi = Mol[i]->phi + k * M;
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (norm);
					if (debug)
					{
						Real sum = lat->ComputeTheta(phi);
						NAMICS_DBG("Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl);
					}
				k++;
			}
		}
	}
	if (charged && neutralizer > -1)
	{
		//}

		//Mol[neutralizer]->phibulk = -A/Mol[neutralizer]->Charge();
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

		if (Mol[solvent]->MolType==water) {

			Mol[solvent]->GetPhib1();
			Mol[solvent]->ComputePhi();
			Mol[solvent]->chainlength=1;
			Mol[solvent]->n = lat->ComputeGN(Mol[solvent]->phi,Mol[solvent]->Markov,M);
			Mol[solvent]->theta=Mol[solvent]->n;
			Mol[solvent]->norm=1.0;

		} else {

			Real norm = Mol[solvent]->phibulk / Mol[solvent]->chainlength;
			Mol[solvent]->n = norm * Mol[solvent]->GN;
			Mol[solvent]->theta = Mol[solvent]->n * Mol[solvent]->chainlength;
			Mol[solvent]->norm = norm;

			int k = 0;

			length = Mol[solvent]->MolMonList.size();
			while (k < length) {
				Real *phi = Mol[solvent]->phi + k * M;
				if (norm > 0)
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (norm);
					if (debug)
					{
						Real sum;
						(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (phi)[__i];
						NAMICS_DBG("Sumphi in mol " << solvent << "for mon " << k << ":" << sum << endl);
					}
				k++;
			}
		}
	}
	if (charged && neutralizer > -1)
	{
		int k = 0;
		length = Mol[neutralizer]->MolMonList.size();
		while (k < length)
		{
			Real *phi = Mol[neutralizer]->phi + k * M;
			if (Mol[neutralizer]->norm > 0)
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] *= (Mol[neutralizer]->norm);
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
		length = SysTagList.size();
		k = 0;
		while (k < length)
		{
			std::copy_n(Seg[SysTagList[k]]->MASK, M, Seg[SysTagList[k]]->phi);
			k++;
		}
		length = SysClampList.size();
		k = 0;
		while (k < length)
		{
			std::copy_n(Seg[SysClampList[k]]->MASK, M, Seg[SysClampList[k]]->phi);
			k++;
		}
	}

	int n_seg = In->MonList.size();
	if (do_blocks) {
		for (int k=0; k<n_mol; k++) {
			if (Mol[k]->freedom =="restricted") {
				Mol[k]->NormPerBlock(split);
				cout<<"*";
			}
		}
		for (int k=0; k<n_seg; k++) {
			if (Seg[k]->freedom =="free") std::fill_n(Seg[k]->phi, M, 0);
		}
		for (int k=0; k<n_mol; k++) {
			int length = Mol[k]->MolMonList.size();
			for (int i=0; i<length; i++) {
				Real *phi_mon = Seg[Mol[k]->MolMonList[i]]->phi;
				Real *phi_molmon = Mol[k]->phi + i * M;
				for (int __i = 0; __i < (M); ++__i) (phi_mon)[__i] += (phi_molmon)[__i];
			}
		}
	}


	for (int i = 0; i < n_seg; i++) {
		lat->set_bounds(Seg[i]->phi);
	}
//(result) = 0; for (int __i = 0; __i < (M); ++__i) (result) += (phitot)[__i];

	if (CalculationType=="steady_state") {
		Real PhiTot0=0;
		Real PhiTotM=0;
		Real Qtot0=0;
		Real QtotM=0;
		int MX=lat->MX;
		int MY=lat->MY;
		int JX=lat->JX;
		int gradients=lat->gradients;

		for (int i = 0; i < n_seg; i++) Seg[i]->PutContraintBC();
					//make sure that both bounds have sumphi=1 and are neutral.

		switch (gradients) {
			case 1:
				for (int i = 0; i < n_seg; i++) {
					if (!(Seg[i]->used_in_mol_nr==solvent || Seg[i]->used_in_mol_nr==neutralizer)) {
						PhiTot0+=Seg[i]->phi[0];
						Qtot0+=Seg[i]->phi[0]*Seg[i]->valence;
						PhiTotM+=Seg[i]->phi[M-1];
						QtotM+=Seg[i]->phi[M-1]*Seg[i]->valence;
					}
				}

				if (Qtot0!=0 && neutralizer==-1) cout <<"Error: neutralizer needed, but was not found. Outcome uncertain" << endl;
				if (Qtot0!=0) {
					Mol[neutralizer]->phitot[0]=-Qtot0/Mol[neutralizer]->Charge(); PhiTot0 +=Mol[neutralizer]->phitot[0];
					Mol[neutralizer]->phitot[M-1]=-QtotM/Mol[neutralizer]->Charge(); PhiTotM +=Mol[neutralizer]->phitot[M-1];
				}

				Mol[solvent]->phitot[0]=1.0-PhiTot0;
				Mol[solvent]->phitot[M-1]=1.0-PhiTotM;

				for (int i=0; i<n_seg; i++) {

					if (Seg[i]->used_in_mol_nr==solvent) {
						Seg[i]->phi[0]=Mol[solvent]->fraction(i)*Mol[solvent]->phitot[0];
						Seg[i]->phi[M-1]=Mol[solvent]->fraction(i)*Mol[solvent]->phitot[M-1];
					}
					if (Seg[i]->used_in_mol_nr==neutralizer) {
						Seg[i]->phi[0]=Mol[neutralizer]->fraction(i)*Mol[neutralizer]->phitot[0];
						Seg[i]->phi[M-1]=Mol[neutralizer]->fraction(i)*Mol[neutralizer]->phitot[M-1];
					}

				}
				break;
			case 2:

					for (int x=1; x<MX+1; x++) {
						PhiTot0=0;
						PhiTotM=0;
						Qtot0=0;
						QtotM=0;
						for (int i = 0; i < n_seg; i++) {
							if (!(Seg[i]->used_in_mol_nr==solvent || Seg[i]->used_in_mol_nr==neutralizer)) {
								PhiTot0+=Seg[i]->phi[x*JX+0];
								Qtot0+=Seg[i]->phi[x*JX+0]*Seg[i]->valence;
								PhiTotM+=Seg[i]->phi[x*JX+MY+1];
								QtotM+=Seg[i]->phi[x*JX+MY+1]*Seg[i]->valence;
							}
						}

						if (Qtot0!=0) {
							Mol[neutralizer]->phitot[x*JX+0]=-Qtot0/Mol[neutralizer]->Charge(); PhiTot0 +=Mol[neutralizer]->phitot[x*JX+0];
							Mol[neutralizer]->phitot[x*JX+MY+1]=-QtotM/Mol[neutralizer]->Charge(); PhiTotM +=Mol[neutralizer]->phitot[x*JX+MY+1];
						}

						Mol[solvent]->phitot[x*JX+0]=1.0-PhiTot0;
						Mol[solvent]->phitot[x*JX+MY+1]=1.0-PhiTotM;

						for (int i=0; i<n_seg; i++) {

							if (Seg[i]->used_in_mol_nr==solvent) {
								Seg[i]->phi[x*JX+0]=Mol[solvent]->fraction(i)*Mol[solvent]->phitot[x*JX+0];
								Seg[i]->phi[x*JX+MY+1]=Mol[solvent]->fraction(i)*Mol[solvent]->phitot[x*JX+MY+1];
							}
							if (Seg[i]->used_in_mol_nr==neutralizer) {
								Seg[i]->phi[x*JX+0]=Mol[neutralizer]->fraction(i)*Mol[neutralizer]->phitot[x*JX+0];
								Seg[i]->phi[x*JX+MY+1]=Mol[neutralizer]->fraction(i)*Mol[neutralizer]->phitot[x*JX+MY+1];
							}

						}
					}
					for (int y=1; y<MY+1; y++) {
						PhiTot0=0;
						PhiTotM=0;
						Qtot0=0;
						QtotM=0;
						for (int i = 0; i < n_seg; i++) {
							if (!(Seg[i]->used_in_mol_nr==solvent || Seg[i]->used_in_mol_nr==neutralizer)) {
								PhiTot0+=Seg[i]->phi[0+y];
								Qtot0+=Seg[i]->phi[0+y]*Seg[i]->valence;
								PhiTotM+=Seg[i]->phi[(MX+1)*JX+y];
								QtotM+=Seg[i]->phi[(MX+1)*JX+y]*Seg[i]->valence;
							}
						}

						if (Qtot0!=0) {
							Mol[neutralizer]->phitot[0+y]=-Qtot0/Mol[neutralizer]->Charge(); PhiTot0 +=Mol[neutralizer]->phitot[0+y];
							Mol[neutralizer]->phitot[(MX+1)*JX+y]=-QtotM/Mol[neutralizer]->Charge(); PhiTotM +=Mol[neutralizer]->phitot[(MX+1)*JX+y];
						}

						Mol[solvent]->phitot[0+y]=1.0-PhiTot0;
						Mol[solvent]->phitot[(MX+1)*JX+y]=1.0-PhiTotM;

						for (int i=0; i<n_seg; i++) {

							if (Seg[i]->used_in_mol_nr==solvent) {
								Seg[i]->phi[0+y]=Mol[solvent]->fraction(i)*Mol[solvent]->phitot[0+y];
								Seg[i]->phi[(MX+1)*JX+y]=Mol[solvent]->fraction(i)*Mol[solvent]->phitot[(MX+1)*JX+y];
							} //else
							if (Seg[i]->used_in_mol_nr==neutralizer) {
								Seg[i]->phi[0+y]=Mol[neutralizer]->fraction(i)*Mol[neutralizer]->phitot[0+y];
								Seg[i]->phi[(MX+1)*JX+y]=Mol[neutralizer]->fraction(i)*Mol[neutralizer]->phitot[(MX+1)*JX+y];
							}

						}
					}
				break;
			default:
				break;
		}


		int it_mon_length=ItMonList.size();
		std::fill_n(B_phitot, M, 0);
		for (int i=0; i<it_mon_length; i++)
			for (int __i = 0; __i < (M); ++__i) (B_phitot)[__i] += (Seg[ItMonList[i]]->B) * (Seg[ItMonList[i]]->phi)[__i];
	}

	for (int i = 0; i < n_seg; i++) {
		Seg[i]->SetPhiSide();
	}


	if (prepare_for_blocks) {
		do_blocks=true;

		for (int k=0; k<n_mol; k++) {
			if (Mol[k]->freedom =="restricted") {
				Mol[k]->SetThetaBlocks(split);
			}
		}
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
	if (CalculationType=="steady_state") {
		CreateMu(lat->M-2); //assuming 1 gradient systems....
		for (int i=0; i<n_mol; i++) {
			Mol[i]->Delta_MU=Mol[i]->Mu;
		}
		CreateMu(1);
		for (int i=0; i<n_mol; i++) {
			Mol[i]->Delta_MU-=Mol[i]->Mu;

		}

	}
	CreateMu(lat->M);

	if ((e_info&& first_pass && CalculationType!="steady_state"))
	{
		cout << "free energy                 = " << FreeEnergy << endl;
		cout << "grand potential             = " << GrandPotential << endl;
	} else {
		if (CalculationType=="steady_state")
		cout << "free energy                 = " << FreeEnergy << endl;
	}
	Real n_times_mu = 0;
	for (int i = 0; i < n_mol; i++)
	{
		Real Mu = Mol[i]->Mu;
		Real n = Mol[i]->n;
		if (Mol[i]->IsClamped())
			n = Mol[i]->n_box;
		n_times_mu += n * Mu;
	}
	if ((e_info && first_pass && CalculationType!="steady_state"))
	{
		cout << "free energy     (GP + n*mu) = " << GrandPotential + n_times_mu << endl;
		cout << "grand potential (F - n*mu)  = " << FreeEnergy - n_times_mu << endl<<endl;;
		//	}

		//	Mol[i]->compute_phi_alias=true;
		//Mol[i]->ComputePhi();
		//}
		//ComputePhis();
		//}
		int M = lat->M;
		for (int i = 0; i < n_mol; i++)
		{
			int n_molmon = Mol[i]->MolMonList.size();
			Real theta_tot = Mol[i]->n * Mol[i]->chainlength;
			for (int j = 0; j < n_molmon; j++)
			{
				Real FRACTION = Mol[i]->fraction(Mol[i]->MolMonList[j]);
				if (Seg[Mol[i]->MolMonList[j]]->freedom != "clamp")
				{
					Real THETA = lat->WeightedSum(Mol[i]->phi + j * M);
					cout << "MOL " << Mol[i]->name << " Fraction " << Seg[Mol[i]->MolMonList[j]]->name << ": " << FRACTION << "=?=" << THETA / theta_tot << " or " << THETA << " of " << theta_tot << endl;
				}
			}
		}
	}
	first_pass=false;

	if (GetValue("compute_Gibbs_excess").size()>0){
		Real Rgibbs=Mol[solvent]->ComputeGibbs(0);
		Real excess=0;
		for (int i=0; i<n_mol; i++) {
			if (i !=solvent) excess+=Mol[i]->ComputeGibbs(Rgibbs);
		}
		if (excess<-1e-5 || excess > 1e-5) cout <<"total excess is not close to zero: " << excess << endl;
	}

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
		//(E) = 0; for (int __i = 0; __i < (M); ++__i) (E) += (phi)[__i] * (side)[__i]; //also need L in nonplaner geometries; for this we need new vector and for the time being I did not do this. (need this result only for one-gradient planar..
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
		if (Mol[i]->freedom == "clamped")
		{
			int n_box = Mol[i]->n_box;
			for (int p=0; p<n_box; p++)
				FreeEnergy -= log(Mol[i]->gn[p]);
		}
		else
		{
			Real n = Mol[i]->n;
			Real GN = Mol[i]->GN;
			int N = Mol[i]->chainlength;
			Real *phi = Mol[i]->phitot; //contains also the tagged segment
			if (Mol[i]->IsTagged())
				N--; //assuming there is just one tagged segment per molecule
			constant = log(N * n / GN) / N;
			std::copy_n(phi, M, TEMP);
			for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (constant);
			for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
		}
	}
	if (Mol[solvent]->MolType==water) Mol[solvent]->AddToF(F);

	if (constraintfields) {
		for (int i=0; i<M; i++) if (beta[i]>0) {
			F[i] +=log(BETA[i])*(Mol[DeltaMolList[0]]->phitot[i]-Mol[DeltaMolList[1]]->phitot[i]);
		}
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
					if (Seg[k]->freedom == "frozen" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "tagged")
						chi = Seg[j]->chi[k];
					else
						chi = Seg[j]->chi[k] / 2; //double counted.
//}

					phi_side = Seg[k]->phi_side;
					if (!(Seg[j]->freedom == "frozen" || Seg[j]->freedom == "clamp" || Seg[j]->freedom == "tagged" || chi == 0))
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
//}

				phi_side = Seg[Sta[i]->mon_nr]->phi_side + Sta[i]->state_nr * M;
				if (!(Seg[j]->freedom == "frozen" || Seg[j]->freedom == "clamp" || Seg[j]->freedom == "tagged" || chi == 0))
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
				if (!(Seg[k]->freedom == "frozen" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "tagged")) chi = chi / 2;
//}


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
//}

			phi_side = Seg[Sta[j]->mon_nr]->phi_side + Sta[j]->state_nr * M;
			if (chi != 0)
			{
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (phi)[__i] * (phi_side)[__i];
				for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (chi);
				for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
			}
		}
		//FL
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
				if (Mol[i]->IsTagged())
				{
					int N = Mol[i]->chainlength;
					if (N > 1)
					{
						fA = fA * N / (N - 1);
						fB = fB * N / (N - 1);
					}
					else
					{
						fA = 0;
						fB = 0;
					}
				}
				Real chi = CHI[Mol[i]->MolMonList[j] * n_mon + Mol[i]->MolMonList[k]] / 2;
				constant -= fA * fB * chi;
			}
		Real *phi = Mol[i]->phitot;
		std::copy_n(phi, M, TEMP);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (constant);
		for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
	}
	//lat->remove_bounds(F);
	for (int __i = 0; __i < (M); ++__i) (F)[__i] = (F)[__i] * (KSAM)[__i]; //clean up contributions in frozen and tagged sites.

std::fill_n(TEMP, M, 0);
	if (charged) {
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] += (q)[__i] * (psi)[__i];
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (0.5);
		for (int __i = 0; __i < (M); ++__i) (F)[__i] += (TEMP)[__i];
	}
	return FreeEnergy + lat->WeightedSum(F);
}

Real System::GetSpontaneousCurvature()
{
	int M = lat->M;
	if (lat->gradients ==1 && lat->geometry=="planar")  return -1.0*lat->MomentPlanar(GrandPotentialDensity,1,M/2+0.5);
	else return 0;
};

Real System::GetKBar()
{
	int M = lat->M;
	if (lat->gradients ==1 && lat->geometry=="planar")  return lat->MomentPlanar(GrandPotentialDensity,2,M/2+0.5);
	else return 0;
};

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
		if (Mol[i]->IsTagged())
		{
			N--;
			phibulk = 0;
		} //One segment of the tagged molecule is tagged and then removed from GP through KSAM
		if (Mol[i]->IsClamped())
		{
			N = N - 2;
			phibulk = 0;
		}
		std::copy_n(phi, M, TEMP);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] = (TEMP)[__i] + (-phibulk);
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (1.0 / N); //GP has wrong sign. will be corrected at end of this routine;
		for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
	}

	if (Mol[solvent]->MolType==water) {Mol[solvent]->AddToGP(GP); }

	for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (alpha)[__i];


	if (constraintfields) {
		for (int i=0; i<M; i++) if (beta[i]>0) {
			GP[i] -=log(BETA[i])*(Mol[DeltaMolList[0]]->phitot[i]-Mol[DeltaMolList[1]]->phitot[i]);
		}
 	}

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
		if (!(Seg[j]->freedom == "tagged" || Seg[j]->freedom == "clamp" || Seg[j]->freedom == "frozen"))
		{
			if (Seg[j]->ns < 2)
			{
				phi = Seg[j]->phi;
				phibulkA = Seg[j]->phibulk;
				for (int k = 0; k < n_seg; k++)
				{
					if (!(Seg[k]->freedom == "tagged" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "frozen"))
					{
						if (Seg[k]->ns < 2)
						{
							chi = Seg[j]->chi[k] / 2;
//}

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
//}

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
			if (!(Seg[k]->freedom == "frozen" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "tagged"))
			{
				if (Seg[k]->ns < 2)
				{
					chi = Sta[j]->chi[k] / 2;
//}

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
//}

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

//out <<"el to G " << lat->WeightedSum(TEMP) << endl; my guess is that I add nothing here....

		for (int __i = 0; __i < M; ++__i) GP[__i] += TEMP[__i];
		for (int __i = 0; __i < M; ++__i) GP[__i] = GP[__i] * KSAM[__i];
		for (int __i = 0; __i < M; ++__i) TEMP[__i] = q[__i] * KSAM[__i];
		for (int __i = 0; __i < M; ++__i) TEMP[__i] = q[__i] - TEMP[__i];
		for (int __i = 0; __i < M; ++__i) TEMP[__i] = TEMP[__i] * psi[__i];
		for (int __i = 0; __i < (M); ++__i) (TEMP)[__i] *= (0.5);
	for (int __i = 0; __i < (M); ++__i) (GP)[__i] += (TEMP)[__i];
}

	if (!charged) for (int __i = 0; __i < (M); ++__i) (GP)[__i] = (GP)[__i] * (KSAM)[__i]; //necessary to make sure that there are no contribution from solid, tagged or clamped sites in GP.

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

		if (Mol[i]->IsTagged())
			NA = NA - 1;
		if (Mol[i]->IsClamped())
			NA = NA - 2;
		if (Mol[i]->IsClamped())
		{
			n = Mol[i]->n_box;
			int n_box = Mol[i]->n_box;
			GN=0;
			for (int p=0; p<n_box; p++)
				GN += log(Mol[i]->gn[p]);
			Mu = -GN / n +1;
		}
		else
		{
			n = Mol[i]->n;
			GN = Mol[i]->GN;
			if (pos==M) Mu = log(NA * n / GN) + 1; else {
				Mu=log(Mol[i]->phitot[pos]) +1;
			}
		}
		constant = 0;
		for (int k = 0; k < n_mol; k++)
		{
			Real NB = Mol[k]->chainlength;
			if (Mol[k]->IsTagged())
				NB = NB - 1;
			if (Mol[k]->IsClamped())
				NB = NB - 2;
			Real phibulkB;
		        if (pos==M) phibulkB=Mol[k]->phibulk; else phibulkB=Mol[k]->phitot[pos];
			constant += phibulkB / NB;
		}
		Mu = Mu - NA * constant;
		if (Mol[solvent]->MolType==water && pos !=M) cout <<"for Moltype==water chemical potential evaluation must be checked in steady state" << endl;
		if (Mol[solvent]->MolType==water) Mu+= NA*(Mol[solvent]->phib1/(1-Mol[solvent]->Kw*Mol[solvent]->phib1)-Mol[solvent]->phibulk);
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
				if (Seg[j]->freedom == "tagged"|| Seg[j]->freedom=="clamped") FA=0; else FA = Mol[i]->fraction(j);
				if (Mol[i]->IsTagged())
					FA *= (NA + 1) / (NA); //works only in case of homopolymers?
				if (Mol[i]->IsClamped())
					FA *= (NA + 2) / (NA); //works only when homopolymers are clamped....needs probably a fix.
				for (int k = 0; k < n_mon; k++)
				{
					if (Seg[k]->ns < 2)
					{
						if (pos==M) phibulkB = Seg[k]->phibulk; else phibulkB=Seg[k]->phi[pos];
						FB = Mol[i]->fraction(k);
						if (Seg[k]->freedom=="tagged"||Seg[k]->freedom=="clamped") FB=0; else FB=Mol[i]->fraction(k);
						if (Mol[i]->IsTagged())
							FB *= (NA + 1) / (NA);
						if (Mol[i]->IsClamped())
							FB *= (NA + 2) / (NA);
						chi = Seg[j]->chi[k] / 2;
						Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
					}
				}
				for (int l = 0; l < statelistlength; l++)
				{     //adjust for steady_state
					phibulkB = Seg[Sta[l]->mon_nr]->state_phibulk[Sta[l]->state_nr];
					FB = Mol[i]->fraction(Sta[l]->mon_nr) * Seg[Sta[l]->mon_nr]->state_alphabulk[Sta[l]->state_nr];
					if (Mol[i]->IsTagged())
						FB *= (NA + 1) / (NA);
					if (Mol[i]->IsClamped())
						FB *= (NA + 2) / (NA);
					chi = Seg[j]->chi[n_mon + l] / 2;
						Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);

				}
			}
		}
		for (int j = 0; j < statelistlength; j++)
		{ //adjust for steady state
			phibulkA = Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
			FA = Mol[i]->fraction(Sta[j]->mon_nr) * Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr];
			if (Mol[i]->IsTagged())
				FA *= (NA + 1) / (NA);
			if (Mol[i]->IsClamped())
				FA *= (NA + 2) / (NA);
			for (int k = 0; k < n_mon; k++)
				if (Seg[k]->ns < 2)
				{
					phibulkB = Seg[k]->phibulk;
					FB = Mol[i]->fraction(k);
					if (Mol[i]->IsTagged())
						FB *= (NA + 1) / (NA);
					if (Mol[i]->IsClamped())
						FB *= (NA + 2) / (NA);
					chi = Sta[j]->chi[k] / 2;
					Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
				}
			for (int k = 0; k < statelistlength; k++)
			{
				phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
				FB = Mol[i]->fraction(Sta[j]->mon_nr) * Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr];
				if (Mol[i]->IsTagged())
					FB *= (NA + 1) / (NA);
				if (Mol[i]->IsClamped())
					FB *= (NA + 2) / (NA);
				chi = Sta[j]->chi[n_mon + k] / 2;
//}

				Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
			}
		}
 

		Mol[i]->Mu = Mu;
	}
	return success;
}
