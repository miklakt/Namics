#include "system.h"
#include "tools.h"
#include <algorithm>
#include <cmath>

System::System(const Input* In_, Lattice* Lat_, std::span<const std::unique_ptr<Segment>> Seg_, std::span<const std::unique_ptr<State>> Sta_, std::span<const std::unique_ptr<Reaction>> Rea_, std::span<const std::unique_ptr<Molecule>> Mol_, std::string name_)
{
	Seg = Seg_;
	Mol = Mol_;
	In = In_;
	name = name_;
	Sta = Sta_;
	Rea = Rea_;
	lat=Lat_;
	NAMICS_DBG( "Constructor for system " << std::endl);
	charged=false;
	grad_epsilon = false;
	all_system=false;
	first_pass=true;
	neutralizer=-1;
}
System::~System()
{
	NAMICS_DBG( "Destructor for system " << std::endl);
	DeAllocateMemory();
}
void System:: DeAllocateMemory(void){
	NAMICS_DBG( "DeAllocateMemory in system " << std::endl);
	if (!all_system) return;
	GrandPotentialDensity.clear();
	FreeEnergyDensity.clear();
	alpha.clear();
	q.clear();
	psi.clear();
	phitot.clear();
	TEMP.clear();
	KSAM.clear();
	CHI.clear();
	EE.clear();
	E.clear();
	psiMask.clear();
	eps.clear();
	all_system=false;
}

void System::AllocateMemory()
{
	NAMICS_DBG( "AllocateMemory in system " << std::endl);
	DeAllocateMemory();
	int M = lat->M;
	GrandPotentialDensity.assign(M, 0);
	FreeEnergyDensity.assign(M, 0);
	alpha.assign(M, 0);
	if (charged)
	{
		q.assign(M, 0);
		psi.assign(M, 0);
		eps.assign(M, 0);
		EE.assign(M, 0);
		E.assign(M, 0);
		psiMask.assign(M, 0);
	}

	phitot.assign(M, 0);
	KSAM.assign(M, 0);
	TEMP.assign(M, 0);
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
	NAMICS_DBG( "generate_mask in system " << std::endl);
	int M = lat->M;
	auto& ksam = KSAM;
	bool success = true;
	FrozenList.clear();
	int length = In->MonList.size();
	for (int i = 0; i < length; i++)
	{

		if (Seg[i]->freedom == "frozen") {
			FrozenList.push_back(i);
		}
	}

	std::fill(ksam.begin(), ksam.end(), 0.0);

	length = FrozenList.size();
	for (int i = 0; i < length; ++i)
	{
		for (int __i = 0; __i < M; ++__i) ksam[__i] += Seg[FrozenList[i]]->MASK[__i];
	}

	for (int __i = 0; __i < M; ++__i) ksam[__i] = ksam[__i] == 0 ? 1 : 0;

	Real accessible_volume = 0;
	if (lat->gradients < 3)
	{
		for (int i = 0; i < M; i++)
		{
			accessible_volume += ksam[i] * lat->L[i];
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
	NAMICS_DBG( "PrepareForCalculations in System " << std::endl);

	bool success = true;
	int M = lat->M;
	auto& ksam = KSAM;
	auto& psi_mask = psiMask;
		success = generate_mask();

	n_mol = In->MolList.size();
	success = lat->PrepareForCalculations();
	int n_mon = In->MonList.size();

	for (int i = 0; i < n_mon; i++)
	{
		success = Seg[i]->PrepareForCalculations(ksam,first_time);
	}
	for (int i = 0; i < n_mol; i++)
	{
		success = Mol[i]->PrepareForCalculations(ksam);
	}
	if (first_time) {
  		if (charged) {
    			int length = FrozenList.size();
    			std::fill(psi_mask.begin(), psi_mask.end(), 0.0);
    			fixedPsi0 = false;
    			for (int i = 0; i < length; ++i) {
      				if (Seg[FrozenList[i]]->fixedPsi0) {
        				fixedPsi0 = true;
        				for (int __i = 0; __i < M; ++__i) psi_mask[__i] += Seg[FrozenList[i]]->MASK[__i];
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
				if (!ContainsValue(SysMonList, Mol[i]->MolMonList[j]))
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
	NAMICS_DBG( "CheckInput for system " << std::endl);
	start=start_;
	bool success = true;
	bool solvent_found = false;
	solvent = -1; //value -1 means no solvent defined.
	Real phibulktot = 0;
	const auto& parameters = In->Parameters("sys", name, start);
	static const std::vector<std::string> keys = {"initial_guess", "guess_inputfile", "write_initial_guess", "X", "E"};
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		success = false;
		std::cout << "sys property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : keys) std::cout << item << std::endl;
	}
	if (success)
	{
		try {
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
				phibulktot=1; std::cout <<"WARNING: no solvent found. Expecting solvent free 'brush'" << std::endl;
			} else {
				std::cout <<"Error: No solvent molecule found. One of your molecules must have 'freedom' : 'solvent'! " << std::endl;
				success=false;
			}
		}
		else
		if (!solvent_found || (phibulktot > 0.99999999 && phibulktot < 1.0000000001))
		{
			std::cout << "In system '" + name + "' the 'solvent' was not found, while the volume fractions of the bulk do not add up to unity. " << std::endl;
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

				if (neutralizer_needed || std::abs(phibulk) > 1e-4)
				{
					std::cout << "Neutralizer needed because some un-pinned molecules have freedom 'restricted' and/or overall charge density of 'free' molecules in the bulk is not equal to zero " << std::endl;
					success = false;
				}
			}
		}
			initial_guess = "previous_result";
			guess_inputfile.clear();
			if (parameters.contains("initial_guess"))
			{
				initial_guess = parameters.at("initial_guess").get<std::string>();
				if (initial_guess != "previous_result" && initial_guess != "file" && initial_guess != "none") {
					std::cout << " Info about 'initial_guess' rejected;" << std::endl;
					success = false;
				}
			}
			if (initial_guess == "file") {
				if (parameters.contains("guess_inputfile")) guess_inputfile = In->ResolvePath(parameters.at("guess_inputfile").get<std::string>());
				else if (In->Start(start).contains("initial_guess")) guess_inputfile = In->json_path;
				else {
					success = false;
					std::cout << "When 'initial_guess' is set to 'file', provide either 'guess_inputfile' or an embedded 'initial_guess' object." << std::endl;
				}
			} else if (initial_guess == "previous_result" && In->Start(start).contains("initial_guess")) {
				initial_guess = "file";
				guess_inputfile = In->json_path;
			}
		write_initial_guess = parameters.value("write_initial_guess", false);
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in system '" << name << "': " << error.what() << std::endl;
			success = false;
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
			std::cout << " num_of_Seg_with_states+num_of_Eqns+num_of_alphabulk_fixed !=num_of_states" << std::endl;
			if (num_of_alphabulk_fixed == 0)
			{
				std::cout << " Consider to define for one of the states an alphabulk value " << std::endl;
			}
			else
			{
				if (num_of_alphabulk_fixed > 1)
				{
					std::cout << " possibly you have specified too many alphabulk values for multiple states " << std::endl;
				}
				else
				{
					if ((num_of_Seg_with_states + num_of_Eqns + num_of_alphabulk_fixed > num_of_states))
					{
						std::cout << " Possibly you have defined too many equations ... " << std::endl;
					}
					else
					{
						std::cout << " Possibly you have defined too few equations ... " << std::endl;
					}
				}
			}
			success = false;
		}
	}

	if (parameters.contains("E"))
	{
		const std::string energy_name = parameters.at("E").get<std::string>();
		if ( !(energy_name=="chi" || energy_name=="Chi" || energy_name=="CHI") ) {
			std::cout <<" Only the FH chi-interactions are implemented. Use 'sys : sysname : E : chi'" << std::endl;
			std::cout <<" Only chi-contributions to E are generated." << std::endl;
		}
	}
	if (parameters.contains("X"))
	{
		XmolList.clear();
		XstateList_1.clear();
		XstateList_2.clear();
		Xn_1.clear();
		std::string s = parameters.at("X").get<std::string>();
		std::vector<std::string> sub;
		std::vector<std::string> SUB;
		In->split(s, '-', sub);
		if (sub[0] != "F")
		{
			std::cout << "X is the characteristic function specified by user." << std::endl;
			std::cout << "Example of how a characteristic function is defined:" << std::endl;
			std::cout << "Case 1: no internal states: 'F - molname_1  -molname_2 - ...' " << std::endl;
			std::cout << "Here molname_1 etc are names of molecules in the system. " << std::endl;
			std::cout << "Case 2: internal states : F - molname_1 -(statename_1,statename_2,#number) - ... " << std::endl;
			std::cout << "Note that in ... you can add as many molname's and (...,...,...) combinations as you wish" << std::endl;
			std::cout << "Here F = Helmholtz energy. " << std::endl;
			std::cout << "and the statement '-molname_i' implies that 'n_i times mu_i' is subtracted from F. " << std::endl;
			std::cout << "The (statename_i,statename_j,n) implies that 'n times theta_i times mu_j' is subtracted from F " << std::endl;
			std::cout << std::endl;
			std::cout << "Error found: The first item is not the expected 'F' " << std::endl;
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
						std::cout << "In characteristic function X, the entry '" + SUB[0] + "' is not a molecule name " << std::endl;
					}
				}
				else
				{ //want to see (statename1,statename2)
					if (SUB.size() != 3)
					{
						std::cout << "In characteristic function X, the entry '" + sub[i] + "' is not recognised as '(statename_1, statename_2,n_1 )' " << std::endl;
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
							std::cout << "In characteristic function X, the entry '" + sub[i] + "' does not include a positive integer at the third argument. " << std::endl;
						}
						else
							Xn_1.push_back(sto);
						if (!(found_1 && found_2))
						{
							std::cout << "In characteristic function X, the entry '" + sub[i] + "' is not coding for '(statename_1,statename_2,n_1)" << std::endl;
							if (!found_1)
								std::cout << "first state name " + SUB[0].substr(1, SUB[0].length() - 1) + " does not exist" << std::endl;
							if (!found_2)
								std::cout << "second state name " + SUB[1] + " does not exist" << std::endl;
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
					std::cout <<"Currently it is not allowed to put a charged frozen segment in boundary. Put this frozen segment inside the system instead. " << std::endl;  success=false;
				}
				bc[Seg[i]->frozen_at_bound]++;
			}
		}
	}
	if (lat->BC[0]=="surface" && bc[0] ==0) {std::cout <<"Lonely 'surface'. Specify a segment with frozen_range including the lowerboundary in x" << std::endl; success=false;}
	if (lat->BC[1]=="surface" && bc[1] ==0) {std::cout <<"Lonely 'surface'. Specify a segment with frozen_range including the lowerboundary in y" << std::endl; success=false;}
	if (lat->BC[2]=="surface" && bc[2] ==0) {std::cout <<"Lonely 'surface'. Specify a segment with frozen_range including the lowerboundary in z" << std::endl; success=false;}
	if (lat->BC[3]=="surface" && bc[3] ==0) {std::cout <<"Lonely 'surface'. Specify a segment with frozen_range including the upperboundary in x" << std::endl; success=false;}
	if (lat->BC[4]=="surface" && bc[4] ==0) {std::cout <<"Lonely 'surface'. Specify a segment with frozen_range including the upperboundary in y" << std::endl; success=false;}
	if (lat->BC[5]=="surface" && bc[5] ==0) {std::cout <<"Lonely 'surface'. Specify a segment with frozen_range including the upperboundary in z" << std::endl; success=false;}
	if (lat->BC[0]=="surface" && bc[0] >1) {std::cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the lowerboundary in x" << std::endl; success=false;}
	if (lat->BC[1]=="surface" && bc[1] >1) {std::cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the lowerboundary in y" << std::endl; success=false;}
	if (lat->BC[2]=="surface" && bc[2] >1) {std::cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the lowerboundary in z" << std::endl; success=false;}
	if (lat->BC[3]=="surface" && bc[3] >1) {std::cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the upperboundary in x" << std::endl; success=false;}
	if (lat->BC[4]=="surface" && bc[4] >1) {std::cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the upperboundary in y" << std::endl; success=false;}
	if (lat->BC[5]=="surface" && bc[5] >1) {std::cout <<"Overpopulated 'surface'. Specify only one segment with frozen_range including the upperboundary in z" << std::endl; success=false;}

	return success;
}

bool System::IsUnique(int Segnr_, int Statenr_)
{
	NAMICS_DBG( "System::IsUnique: Segnr = " << Segnr_ << " Statenr = " << Statenr_ << std::endl);
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

void System::PushOutput()
{
	NAMICS_DBG( "PushOutput for system " << std::endl);
	const auto& parameters = In->Parameters("sys", name, start);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["e"] = e;
	OUTPUT["k_B"] = k_B;
	OUTPUT["eps0"] = eps0;
	OUTPUT["temperature"] = T;
	OUTPUT["free_energy"] = FreeEnergy;
	OUTPUT["grand_potential"] = GrandPotential;
	OUTPUT["start"] = start;
	if (lat->gradients==1) {
		OUTPUT["Laplace_pressure"] = -GrandPotentialDensity[lat->fjc];
	}
	if (lat->gradients==2) {
		if (lat->BC[4]=="surface") {
			OUTPUT["Laplace_pressure"] = -GrandPotentialDensity[lat->P(2*lat->fjc,(lat->MY+lat->fjc)/2)];
		} else {
			OUTPUT["Laplace_pressure"] = -GrandPotentialDensity[lat->P(2*lat->fjc,lat->MY)];
		}
	}
	if (parameters.contains("E"))
	{
		Real sumE=0;
		int length= In->MonList.size();
		for (int i=0; i<length; i++)
		for (int j=i+1; j<length; j++) {
			Real Eij=GetE(i,j);
			OUTPUT["I_" + Seg[i]->name + "_" + Seg[j]->name] = Eij;
			OUTPUT["I_" + Seg[j]->name + "_" + Seg[i]->name] = Eij;
			sumE+=Eij*Seg[i]->chi[j];
		}
		OUTPUT["E"] = sumE;
	}
	int n_seg=In->MonList.size();
	for (int i=0; i<n_seg; i++)
	for (int j=0; j<n_seg; j++){
		OUTPUT["chi_" + Seg[i]->name + "_" + Seg[j]->name] = CHI[i * n_seg + j];
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
				std::cout << "Failed to find chemical potential for state " + In->StateList[XstateList_2[i]] + ": (not a monomer?) In characteristic function X, mu is set to zero." << std::endl;
				mu = 0;
			}
			X -= Seg[Sta[XstateList_1[i]]->mon_nr]->state_theta[Sta[XstateList_1[i]]->state_nr] * Xn_1[i] * mu;
		}
		OUTPUT["X"] = X;
		std::cout << " X  = " << X << std::endl;
	}
	if (solvent > -1) OUTPUT["solvent"] = Mol[solvent]->name;
	OUTPUT["alpha"] = {{"profile", 0}};
	OUTPUT["GrandPotentialDensity"] = {{"profile", 1}};
	OUTPUT["grand_potential_density"] = {{"profile", 1}};
	OUTPUT["FreeEnergyDensity"] = {{"profile", 2}};
	OUTPUT["free_energy_density"] = {{"profile", 2}};
	OUTPUT["phitot"] = {{"profile", 6}};

	if (charged)
	{
		OUTPUT["Dpsi"] = psi[lat->M-1] - psi[0];
		OUTPUT["psi"] = {{"profile", 3}};
		OUTPUT["q"] = {{"profile", 4}};
		OUTPUT["eps"] = {{"profile", 5}};
	}
}

std::span<Real> System::GetPointer(int profile)
{
	NAMICS_DBG( "GetPointer for system " << std::endl);
	if (profile == 0)
		return alpha;
	if (profile == 1)
		return GrandPotentialDensity;
	if (profile == 2)
		return FreeEnergyDensity;
	if (profile == 3)
		return psi;
	if (profile == 4)
		return q;
	if (profile == 5)
		return eps;
	if (profile == 6)
		return phitot;
	return {};
}

bool System::CheckChi_values(int n_seg)
{
	NAMICS_DBG( "CheckChi_values for system " << std::endl);
	bool success = true;
	const int n_segments = n_seg;
	const int n_states = In->StateList.size();
	const int n_chi = n_segments + n_states;
	CHI.assign(n_segments * n_segments, 0);
	for (int i = 0; i < n_segments; i++) Seg[i]->chi.assign(n_chi, 0);
	for (int i = 0; i < n_states; i++) Sta[i]->chi.assign(n_chi, 0);

	auto read_chi = [&](const ParameterStore& parameters, const std::string& target, Real& value) {
		const std::string chi_key = "chi_" + target;
		if (!parameters.contains(chi_key)) return false;
		value = parameters.at(chi_key).get<Real>();
		return true;
	};

	for (int i = 0; i < n_segments; i++) {
		const auto& left_parameters = In->Parameters("mon", Seg[i]->name, start);
		for (int j = i; j < n_segments; j++) {
			Real left = 0;
			Real right = 0;
			const bool has_left = read_chi(left_parameters, Seg[j]->name, left);
			const bool has_right = read_chi(In->Parameters("mon", Seg[j]->name, start), Seg[i]->name, right);
			Real chi = 0;
			if (i != j) {
				if (has_left && has_right && left != right) {
					success = false;
					std::cout << " conflict in chi values! chi(" << Seg[i]->name << "," << Seg[j]->name << ") != chi (" << Seg[j]->name << "," << Seg[i]->name << ")" << std::endl;
					chi = left;
				} else if (has_left) {
					chi = left;
				} else if (has_right) {
					chi = right;
				}
			}
			CHI[i * n_segments + j] = CHI[j * n_segments + i] = chi;
		}
	}

	std::vector<int> copy_source(n_segments);
	for (int i = 0; i < n_segments; i++) {
		copy_source[i] = i;
		const std::string copy_of = In->Parameters("mon", Seg[i]->name, start).value("set_equal_to", std::string{});
		if (copy_of.empty()) continue;
		int segnr = -1;
		for (int j = 0; j < n_segments; j++) {
			if (Seg[j]->name == copy_of) {
				segnr = j;
				break;
			}
		}
		if (segnr < 0 || segnr == i) {
			if (segnr < 0) std::cout <<"In segment " << Seg[i]->name << " 'set_equal_to' is rejected because the segment " << copy_of << " was not found" << std::endl;
			else std::cout <<"In segment " << Seg[i]->name << " 'set_equal_to' is rejected because the segment " << copy_of << " can not copied from itself...." << std::endl;
			continue;
		}
		copy_source[i] = segnr;
		Seg[i]->epsilon = Seg[segnr]->epsilon;
		for (int j = 0; j < n_segments; j++) {
			CHI[i * n_segments + j] = CHI[segnr * n_segments + j];
			CHI[j * n_segments + i] = CHI[i * n_segments + j];
		}
		CHI[i * n_segments + segnr] = CHI[segnr * n_segments + i] = 0;
		CHI[i * n_segments + i] = 0;
	}

	for (int i = 0; i < n_segments; i++) {
		std::copy_n(CHI.begin() + i * n_segments, n_segments, Seg[i]->chi.begin());
	}

	for (int i = 0; i < n_states; i++) {
		const auto& left_parameters = In->Parameters("state", Sta[i]->name, start);
		for (int j = i; j < n_states; j++) {
			Real left = 0;
			Real right = 0;
			const bool has_left = read_chi(left_parameters, Sta[j]->name, left);
			const bool has_right = read_chi(In->Parameters("state", Sta[j]->name, start), Sta[i]->name, right);
			Real chi = i == j ? 0 : Seg[Sta[i]->mon_nr]->chi[Sta[j]->mon_nr];
			if (i != j) {
				if (has_left && has_right && left != right) {
					success = false;
					std::cout << " conflict in chi values! chi(" << Sta[i]->name << "," << Sta[j]->name << ") != chi (" << Sta[j]->name << "," << Sta[i]->name << ")" << std::endl;
					chi = left;
				} else if (has_left) {
					chi = left;
				} else if (has_right) {
					chi = right;
				}
			}
			Sta[i]->chi[n_segments + j] = Sta[j]->chi[n_segments + i] = chi;
		}
	}

	for (int i = 0; i < n_segments; i++) {
		const auto& left_parameters = In->Parameters("mon", Seg[copy_source[i]]->name, start);
		for (int j = 0; j < n_states; j++) {
			Real left = 0;
			Real right = 0;
			const bool has_left = read_chi(left_parameters, Sta[j]->name, left);
			const bool has_right = read_chi(In->Parameters("state", Sta[j]->name, start), Seg[i]->name, right);
			Real chi = Seg[i]->chi[Sta[j]->mon_nr];
			if (has_left && has_right && left != right) {
				success = false;
				std::cout << " conflict in chi values! chi(" << Seg[i]->name << "," << Sta[j]->name << ") != chi (" << Sta[j]->name << "," << Seg[i]->name << ")" << std::endl;
				chi = left;
			} else if (has_left) {
				chi = left;
			} else if (has_right) {
				chi = right;
			}
			if (Sta[j]->mon_nr == i && chi != 0) {
				success = false;
				std::cout << " chi between mon-type and one of its states is not allowed for chi(" << Sta[j]->name << "," << Seg[i]->name << ")" << std::endl;
				chi = 0;
			}
			Seg[i]->chi[n_segments + j] = Sta[j]->chi[i] = chi;
		}
	}

	return success;
}

void System::DoElectrostatics(std::span<Real> g, std::span<const Real> x)
{
	int M = lat->M;
	int n_seg = In->MonList.size();
	auto& q = this->q;
	auto& eps = this->eps;
	auto& phitot = this->phitot;
	auto& psi = this->psi;
	std::fill(q.begin(), q.end(), 0.0);
	std::fill(eps.begin(), eps.end(), 0.0);
	for (int i = 0; i < n_seg; i++)
	{
		if (Seg[i]->ns < 2)
		{
			for (int __i = 0; __i < M; ++__i) q[__i] += Seg[i]->valence * Seg[i]->phi[__i];
		}
		lat->set_bounds(Seg[i]->phi.data());
		for (int __i = 0; __i < M; ++__i) eps[__i] += Seg[i]->epsilon * Seg[i]->phi[__i];
	}
	int statelistlength = In->StateList.size();
	for (int i = 0; i < statelistlength; i++)
	{
		auto phi_state = std::span<const Real>(Seg[Sta[i]->mon_nr]->phi_state.data() + Sta[i]->state_nr * M,
		                                      static_cast<size_t>(M));
		for (int __i = 0; __i < M; ++__i) q[__i] += Sta[i]->valence * phi_state[__i];
	}

	for (int __i = 0; __i < M; ++__i) q[__i] = phitot[__i] != 0 ? q[__i] / phitot[__i] : 0;
	for (int __i = 0; __i < M; ++__i) eps[__i] = phitot[__i] != 0 ? eps[__i] / phitot[__i] : 0;
	std::copy(x.begin(), x.end(), psi.begin());
	std::copy(psi.begin(), psi.end(), g.begin());
	lat->set_M_bounds(psi.data());
	if (fixedPsi0) {
		int length=FrozenList.size();
		for (int i=0; i<length; i++) {
			Seg[FrozenList[i]]->UpdateValence(g.data(),psi,q,eps,grad_epsilon);
		}
	}
}

void System:: ComputePhis(std::span<const Real> x,bool first_time, Real residual, bool final_pass) {
NAMICS_DBG("ComputPhis in  system " << std::endl);
	PutU(x);
	PrepareForCalculations(first_time);
	ComputePhis(residual, final_pass);
}

bool System:: PutU(std::span<const Real> xx) {
NAMICS_DBG("PutU in  Solve " << std::endl);
	int M=lat->M;
	auto& psi = this->psi;
	auto& ee = EE;
	int itmonlistlength=ItMonList.size();
	int itstatelistlength=ItStateList.size();
	int monlistlength =In->MonList.size();
	int statelistlength=In->StateList.size();
	int k=0;

	int itpos=(itmonlistlength+itstatelistlength)*M;
	Real valence;
	bool success=true;
	const auto add_segment_contributions = [&](std::span<Real> field, const Segment& seg, Real segment_valence) {
		for (int __i = 0; __i < M; ++__i) field[__i] += seg.u_ext[__i];
		if (charged){
			for (int __i = 0; __i < M; ++__i) field[__i] += (-1.0 * seg.epsilon) * ee[__i];
			if (segment_valence !=0) {
				for (int __i = 0; __i < M; ++__i) field[__i] += segment_valence * psi[__i];
			}
		}
	};

	if (charged) {
		auto psi_input = xx.subspan(static_cast<size_t>(itpos), static_cast<size_t>(M));
		std::copy(psi_input.begin(), psi_input.end(), psi.begin());
		lat->UpdateEE(ee.data(),psi.data(),E.data());
	}


	for (int i=0; i<itmonlistlength; i++) {
		int IM=ItMonList[i];
		auto u = std::span<Real>(Seg[IM]->u.data(), static_cast<size_t>(M));
		std::copy(xx.begin() + k*M, xx.begin() + (k+1)*M, u.begin());
		valence=Seg[IM]->valence;
		add_segment_contributions(u, *Seg[IM], valence);
		for (int j=0; j<monlistlength; j++) {
			if (Seg[j]->seg_nr_of_copy==IM && Seg[j]->ns<2) {
				u = std::span<Real>(Seg[j]->u.data(), static_cast<size_t>(M));
				std::copy(xx.begin() + k*M, xx.begin() + (k+1)*M, u.begin());
				valence=Seg[j]->valence;
				add_segment_contributions(u, *Seg[j], valence);
			}

		}
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->seg_nr_of_copy==IM) {
				u = std::span<Real>(Seg[Sta[j]->mon_nr]->u.data()+Sta[j]->state_nr*M, static_cast<size_t>(M));
				std::copy(xx.begin() + k*M, xx.begin() + (k+1)*M, u.begin());
				valence=Sta[j]->valence;
				add_segment_contributions(u, *Seg[Sta[j]->mon_nr], valence);
			}
		}
		k++;
	}

	for (int i=0; i<itstatelistlength; i++) {
		int IS=ItStateList[i];
		auto u = std::span<Real>(Seg[Sta[IS]->mon_nr]->u.data()+(Sta[IS]->state_nr)*M, static_cast<size_t>(M));
		std::copy(xx.begin() + k*M, xx.begin() + (k+1)*M, u.begin());
		valence=Sta[IS]->valence;
		add_segment_contributions(u, *Seg[Sta[IS]->mon_nr], valence);
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->state_nr_of_copy==IS) {
				u = std::span<Real>(Seg[Sta[j]->mon_nr]->u.data()+Sta[j]->state_nr*M, static_cast<size_t>(M));
				std::copy(xx.begin() + k*M, xx.begin() + (k+1)*M, u.begin());
				valence=Sta[j]->valence;
				add_segment_contributions(u, *Seg[Sta[j]->mon_nr], valence);
			}
		}
		k++;
	}
	if (charged) itpos +=M;

	return success;
}

void System::Classical_residual(std::span<const Real> x,std::span<Real> g,Real residual, int iterations){
NAMICS_DBG("Classical_residuals in scf mode in system " << std::endl);
	int M=lat->M;
	auto& phitot = this->phitot;
	auto& alpha = this->alpha;
	auto& ksam = this->KSAM;
	auto& psi = this->psi;
	auto& q = this->q;
	auto& eps = this->eps;
	auto& psi_mask = this->psiMask;
	Real chi;
	int mon_length = In->MonList.size(); //also frozen segments
	int i,k;

	int itmonlistlength=ItMonList.size();
	int state_length = In->StateList.size();
	int itstatelistlength=ItStateList.size();

	std::copy(x.begin(), x.end(), g.begin());
	ComputePhis(x,iterations==0,residual,false);
 	std::fill(alpha.begin(), alpha.end(), 0.0);

	for (i=0; i<itmonlistlength; i++) {
		auto g_mon = g.subspan(static_cast<size_t>(i * M), static_cast<size_t>(M));
		for (int __i = 0; __i < M; ++__i) g_mon[__i] += Seg[ItMonList[i]]->u_ext[__i];
		for (k=0; k<mon_length; k++) {
			if (Seg[k]->ns<2) {
				chi =Seg[ItMonList[i]]->chi[k];
				if (chi!=0) {
					for (int __i = 0; __i < M; ++__i) if (phitot[__i] > 0) g_mon[__i] -= chi * ((Seg[k]->phi_side[__i] / phitot[__i]) - Seg[k]->phibulk);
				}
			}
		}
		for (k=0; k<state_length; k++) {
			chi =Seg[ItMonList[i]]->chi[mon_length+k];
			if (chi!=0) {
				auto phi_side = std::span<const Real>(Seg[Sta[k]->mon_nr]->phi_side).subspan(static_cast<size_t>(Sta[k]->state_nr * M), static_cast<size_t>(M));
				for (int __i = 0; __i < M; ++__i) if (phitot[__i] > 0) g_mon[__i] -= chi * ((phi_side[__i] / phitot[__i]) - Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr]);
			}
		}
	}
	for (i=0; i<itmonlistlength; i++) {
		auto g_mon = g.subspan(static_cast<size_t>(i * M), static_cast<size_t>(M));
		for (int __i = 0; __i < M; ++__i) alpha[__i] += g_mon[__i];
	}

	for (i=0; i<itstatelistlength; i++) {
		auto g_state = g.subspan(static_cast<size_t>((itmonlistlength + i) * M), static_cast<size_t>(M));
		for (k=0; k<mon_length; k++) {
			if (Seg[k]->ns<2) {
				chi =Sta[ItStateList[i]]->chi[k];
				if (chi!=0) {
					for (int __i = 0; __i < M; ++__i) if (phitot[__i] > 0) g_state[__i] -= chi * ((Seg[k]->phi_side[__i] / phitot[__i]) - Seg[k]->phibulk);
				}
			}
		}


		for (k=0; k<state_length; k++) {
			chi =Sta[ItStateList[i]]->chi[mon_length+k];
			if (chi!=0) {
				auto phi_side = std::span<const Real>(Seg[Sta[k]->mon_nr]->phi_side).subspan(static_cast<size_t>(Sta[k]->state_nr * M), static_cast<size_t>(M));
				for (int __i = 0; __i < M; ++__i) if (phitot[__i] > 0) g_state[__i] -= chi * ((phi_side[__i] / phitot[__i]) - Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr]);
			}
		}
	}
	for (i=0; i<itstatelistlength; i++) {
		auto g_state = g.subspan(static_cast<size_t>((itmonlistlength + i) * M), static_cast<size_t>(M));
		for (int __i = 0; __i < M; ++__i) alpha[__i] += g_state[__i];
	}
	for (int __i = 0; __i < M; ++__i) alpha[__i] *= (1.0/(itmonlistlength+itstatelistlength));
	for (i=0; i<itmonlistlength; i++) {
		auto g_mon = g.subspan(static_cast<size_t>(i * M), static_cast<size_t>(M));
		for (int __i = 0; __i < M; ++__i) g_mon[__i] = phitot[__i] > 0 ? g_mon[__i] - alpha[__i] + 1 / phitot[__i] - 1.0 : 0;
		lat->remove_bounds(g_mon.data());
		for (int __i = 0; __i < M; ++__i) g_mon[__i] *= ksam[__i];
	}
	for (i=0; i<itstatelistlength; i++) {
		auto g_state = g.subspan(static_cast<size_t>((itmonlistlength + i) * M), static_cast<size_t>(M));
		for (int __i = 0; __i < M; ++__i) g_state[__i] = phitot[__i] > 0 ? g_state[__i] - alpha[__i] + 1 / phitot[__i] - 1.0 : 0;
		lat->remove_bounds(g_state.data());
		for (int __i = 0; __i < M; ++__i) g_state[__i] *= ksam[__i];
	}

	int itpos=(itmonlistlength+itstatelistlength)*M;

	if (charged) {
		auto x_psi = x.subspan(static_cast<size_t>(itpos), static_cast<size_t>(M));
		auto g_psi = g.subspan(static_cast<size_t>(itpos), static_cast<size_t>(M));
		std::copy(x_psi.begin(), x_psi.end(), g_psi.begin());
		DoElectrostatics(g_psi,x_psi);
		lat->set_M_bounds(psi.data());
		if (M > 1) psi[0] = psi[1];
		lat->UpdatePsi(g_psi.data(),psi.data(),q.data(),eps.data(),psi_mask.data(),grad_epsilon,fixedPsi0);
		lat->remove_bounds(g_psi.data());
	}
}


bool System::ComputePhis(Real residual, bool final_pass){
NAMICS_DBG("ComputePhis in system" << std::endl);
	int M= lat->M;
	const auto slice_len = static_cast<size_t>(M);
	auto& phitot = this->phitot;
	const auto normalize_ranked_phi = [&](Molecule& mol, Real norm, bool divide_by_g1) {
		if (!final_pass || mol.phi_ranked.empty()) return;
		int s = 0;
		for (size_t b = 0; b < mol.mon_nr.size(); ++b) {
			auto g1 = std::span<const Real>(Seg[mol.mon_nr[b]]->G1);
			for (int k = 0; k < mol.n_mon[b]; ++k, ++s) {
				auto phi = std::span<Real>(mol.phi_ranked).subspan(static_cast<size_t>(s * M), slice_len);
				if (divide_by_g1) {
					for (int __i = 0; __i < M; ++__i) phi[__i] = g1[__i] != 0 ? phi[__i] / g1[__i] : 0;
				}
				if (norm > 0) {
					for (int __i = 0; __i < M; ++__i) phi[__i] *= norm;
				}
			}
		}
	};
	const auto normalize_molecule = [&](Molecule& mol, Real norm, bool normalize_ranked, bool divide_by_g1) {
		if (mol.freedom == "frozen") return;
		for (size_t k = 0; k < mol.MolMonList.size(); ++k) {
			const int seg = mol.MolMonList[k];
			auto phi = std::span<Real>(mol.phi).subspan(static_cast<size_t>(k * M), slice_len);
			auto g1 = std::span<const Real>(Seg[seg]->G1);
			if (divide_by_g1) {
				for (int __i = 0; __i < M; ++__i) phi[__i] = g1[__i] != 0 ? phi[__i] / g1[__i] : 0;
			}
			if (norm > 0) {
				for (int __i = 0; __i < M; ++__i) phi[__i] *= norm;
			}
		}
		if (normalize_ranked) normalize_ranked_phi(mol, norm, divide_by_g1);
	};
	Real A=0, B=0; //A should contain sum_phi*charge; B should contain sum_phi
	bool success=true;
	std::fill(phitot.begin(), phitot.end(), 0.0);

	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		for (int __i = 0; __i < M; ++__i) phitot[__i] += Seg[FrozenList[i]]->phi[__i];
	}

	for (int i=0; i<n_mol; i++) {
		if (final_pass && !Mol[i]->phi_ranked.empty()) {
			std::fill(Mol[i]->phi_ranked.begin(), Mol[i]->phi_ranked.end(), 0);
		}
		success = Mol[i]->ComputePhi(final_pass);
	}

	for (int i = 0; i < n_mol; i++) {
		auto& mol = *Mol[i];
		Real norm = 0;
		if (mol.freedom == "free") {
			norm = mol.phibulk / mol.chainlength;
			mol.n = norm * mol.GN;

			A += mol.phibulk * mol.Charge();
			B += mol.phibulk;
		} else if (mol.freedom == "restricted") {
			if (mol.GN > 0) {
				norm = mol.n / mol.GN;
				if (mol.IsPinned()) {
					mol.phibulk = 0;
				} else {
					mol.phibulk = mol.chainlength * norm;
					A += mol.phibulk * mol.Charge();
					B += mol.phibulk;
				}
			} else {
				norm = 0;
				std::cout << "GN for molecule " << i << " is not larger than zero..." << std::endl;
				std::cout << "Consider compiling with LongReal enabled if this is an overflow issue." << std::endl;
				throw - 1;
			}
		}

		if (mol.IsPinned()) {
			if (mol.GN > 0) {
				if (mol.n == 0) mol.n = 1;
				norm = mol.n / mol.GN;
			} else {
				norm = 0;
				std::cout << "GN for molecule " << i << " is not larger than zero..." << std::endl;
			}
			mol.phibulk = 0;
		}
		mol.norm = norm;
		normalize_molecule(mol, norm, true, true);
	}
	if (charged && neutralizer > -1) {
		auto& neutral = *Mol[neutralizer];
		const Real solvent_charge = Mol[solvent]->Charge();
		const Real neutral_charge = neutral.Charge();
		if (neutral_charge == solvent_charge) {
			std::cout << "WARNING: solvent charge equals neutralizer charge; outcome problematic...." << std::endl;
		} else {
			neutral.phibulk = ((B - 1.0) * solvent_charge - A) / (neutral_charge - solvent_charge);
		}

		if (neutral.phibulk < 0) {
			std::cout << "WARNING: neutralizer has negative phibulk. Consider changing neutralizer...: outcome problematic...." << std::endl;
			std::cout << "A is " << A << std::endl;
			for (int j = 0; j < n_mol; j++) {
				std::cout << " mol : " << Mol[j]->name << " phibulk " << Mol[j]->phibulk << std::endl;
			}
		}

		B += neutral.phibulk;
		Real norm = neutral.phibulk / neutral.chainlength;
		neutral.n = norm * neutral.GN;
		neutral.theta = neutral.n * neutral.chainlength;
		neutral.norm = norm;
		normalize_molecule(neutral, norm, true, false);
	}

	if (solvent > -1) {
		auto& solvent_mol = *Mol[solvent];
		solvent_mol.phibulk = 1.0 - B;
		if (solvent_mol.phibulk < 0) {
			std::cout << "WARNING: solvent has negative phibulk. outcome problematic " << std::endl;
			throw - 4;
		}

		Real norm = solvent_mol.phibulk / solvent_mol.chainlength;
		solvent_mol.n = norm * solvent_mol.GN;
		solvent_mol.theta = solvent_mol.n * solvent_mol.chainlength;
		solvent_mol.norm = norm;
		normalize_molecule(solvent_mol, norm, true, false);
	}

	for (int i = 0; i < n_mol; i++) {
		auto& mol = *Mol[i];
		for (size_t k = 0; k < mol.MolMonList.size(); ++k) {
			const int seg = mol.MolMonList[k];
			auto& phi_mon = Seg[seg]->phi;
			auto& mol_phitot = mol.phitot;
			auto phi_molmon = std::span<const Real>(mol.phi).subspan(static_cast<size_t>(k * M), slice_len);
			for (int __i = 0; __i < M; ++__i) phi_mon[__i] += phi_molmon[__i];
			for (int __i = 0; __i < M; ++__i) phitot[__i] += phi_molmon[__i];
			for (int __i = 0; __i < M; ++__i) mol_phitot[__i] += phi_molmon[__i];
			Seg[seg]->phibulk += mol.fraction(seg) * mol.phibulk;
		}
	}

	int n_seg = In->MonList.size();
	for (int i = 0; i < n_seg; i++) {
		lat->set_bounds(Seg[i]->phi.data());
	}

	for (int i = 0; i < n_seg; i++) {
		Seg[i]->SetPhiSide();
	}
	return success;
}

bool System::CheckResults(bool e_info_)
{
	NAMICS_DBG( "CheckResults for system " << std::endl);

	bool e_info = e_info_;
	bool success = true;

	FreeEnergy = GetFreeEnergy();
	GrandPotential = GetGrandPotential();
	CreateMu(lat->M);
	int M = lat->M;
	for (int i = 0; i < n_mol; ++i) {
		if (!std::isfinite(Mol[i]->n)) {
			std::cerr << "Detected invalid n in computed solver state for molecule " << Mol[i]->name << "." << std::endl;
			return false;
		}
		if (!std::isfinite(Mol[i]->theta)) {
			std::cerr << "Detected invalid theta in computed solver state for molecule " << Mol[i]->name << "." << std::endl;
			return false;
		}
		if (!std::isfinite(Mol[i]->phibulk)) {
			std::cerr << "Detected invalid phibulk in computed solver state for molecule " << Mol[i]->name << "." << std::endl;
			return false;
		}
		if (!std::isfinite(Mol[i]->GN)) {
			std::cerr << "Detected invalid GN in computed solver state for molecule " << Mol[i]->name << "." << std::endl;
			return false;
		}
		for (int j = 0; j < M; ++j) {
			if (!std::isfinite(Mol[i]->phitot[j])) {
				std::cerr << "Detected invalid phitot in computed solver state for molecule " << Mol[i]->name << "." << std::endl;
				return false;
			}
		}
		const int n_molmon = Mol[i]->MolMonList.size();
		for (int j = 0; j < n_molmon * M; ++j) {
			if (!std::isfinite(Mol[i]->phi[j])) {
				std::cerr << "Detected invalid phi in computed solver state for molecule " << Mol[i]->name << "." << std::endl;
				return false;
			}
		}
	}

	if (e_info && first_pass)
	{
		std::cout << "free energy                 = " << FreeEnergy << std::endl;
		std::cout << "grand potential             = " << GrandPotential << std::endl;
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
		std::cout << "free energy     (GP + n*mu) = " << GrandPotential + n_times_mu << std::endl;
		std::cout << "grand potential (F - n*mu)  = " << FreeEnergy - n_times_mu << std::endl<<std::endl;;
		for (int i = 0; i < n_mol; i++)
		{
			int n_molmon = Mol[i]->MolMonList.size();
			Real theta_tot = Mol[i]->n * Mol[i]->chainlength;
			for (int j = 0; j < n_molmon; j++)
			{
				Real FRACTION = Mol[i]->fraction(Mol[i]->MolMonList[j]);
				auto phi = std::span<Real>(Mol[i]->phi).subspan(static_cast<size_t>(j * M), static_cast<size_t>(M));
				Real THETA = lat->WeightedSum(phi.data());
				std::cout << "MOL " << Mol[i]->name << " Fraction " << Seg[Mol[i]->MolMonList[j]]->name << ": " << FRACTION << "=?=" << THETA / theta_tot << " or " << THETA << " of " << theta_tot << std::endl;
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
	std::vector<Real> temp(M, 0);
	const auto& L = lat->L;
	const auto& phi = Seg[Seg1]->phi;
	const auto& side = Seg[Seg2]->phi_side;
	if (Seg1!=Seg2) {
		for (int __i = 0; __i < M; ++__i) temp[__i] = L[__i] * phi[__i];
		for (int __i = 0; __i < M; ++__i) temp[__i] *= side[__i];
		(E) = 0; for (int __i = 0; __i < M; ++__i) (E) += temp[__i];
	}
	return E;//only the contacts
}

Real System::GetFreeEnergy(void)
{ //eqn 2.91 of thesis of J.v.Male;
	NAMICS_DBG( "GetFreeEnergy for system " << std::endl);
	int M = lat->M;
	Real FreeEnergy = 0;
	auto& F = FreeEnergyDensity;
	auto& temp = TEMP;
	auto& q = this->q;
	auto& psi = this->psi;
	auto& ksam = KSAM;
	Real constant = 0;
	int n_seg = In->MonList.size();
	int n_mol = In->MolList.size();
	int n_states = In->StateList.size();
	for (int i=0; i<n_mol; i++) {
		lat->remove_bounds(Mol[i]->phitot.data());
	}
	int n_mon = In->MonList.size();
	for (int i = 0; i < n_mon; i++)
	{
		if (Seg[i]->ns < 2)
			lat->remove_bounds(Seg[i]->phi_side.data());
		else
		{
			for (int j = 0; j < Seg[i]->ns; j++)
			{
				lat->remove_bounds(Seg[i]->phi_side.data() + j * M);
			}
		}
	}

	std::fill(F.begin(), F.end(), 0.0);

	for (int i = 0; i < n_mol; i++)
	{
		Real n = Mol[i]->n;
		Real GN = Mol[i]->GN;
		int N = Mol[i]->chainlength;
		auto phi = std::span<const Real>(Mol[i]->phitot);
		constant = std::log(N * n / GN) / N;
		std::copy(phi.begin(), phi.end(), temp.begin());
		for (int __i = 0; __i < M; ++__i) temp[__i] *= constant;
		for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
	}
	Real chi;
	int n_sysmon = SysMonList.size();
	for (int j = 0; j < n_sysmon; j++)
	{
		int n_states=Seg[SysMonList[j]]->ns;

		if (n_states==1) {
			auto phi = std::span<const Real>(Seg[SysMonList[j]]->phi);
			auto g = std::span<const Real>(Seg[SysMonList[j]]->G1);
			for (int __i = 0; __i < M; ++__i) temp[__i] = g[__i] > 0 ? -std::log(g[__i]) : 0;
			for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * temp[__i];
			for (int __i = 0; __i < M; ++__i) temp[__i] *= -1;
			for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
		} else {
			for (int k=0; k<n_states; k++) {
				auto phi = std::span<const Real>(Seg[SysMonList[j]]->phi_state).subspan(static_cast<size_t>(k * M), static_cast<size_t>(M));
				auto u = std::span<const Real>(Seg[SysMonList[j]]->u).subspan(static_cast<size_t>(k * M), static_cast<size_t>(M));
				std::copy(u.begin(), u.end(), temp.begin());
				for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * temp[__i];
				for (int __i = 0; __i < M; ++__i) temp[__i] *= -1;
				for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
			}
		}
	}

	for (int j = 0; j < n_seg; j++)
	{
		if (Seg[j]->ns < 2)
		{
			auto phi = std::span<const Real>(Seg[j]->phi);
			for (int k = 0; k < n_seg; k++)
			{
				if (Seg[k]->ns < 2)
				{
					if (Seg[k]->freedom == "frozen")
						chi = Seg[j]->chi[k];
					else
						chi = Seg[j]->chi[k] / 2; //double counted.

					auto phi_side = std::span<const Real>(Seg[k]->phi_side);
					if (!(Seg[j]->freedom == "frozen" || chi == 0))
					{
						for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
						for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
						for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
					}
				}
			}
			for (int i = 0; i < n_states; i++)
			{
					chi = Seg[j]->chi[n_seg + i]/2;

					auto phi_side = std::span<const Real>(Seg[Sta[i]->mon_nr]->phi_side).subspan(static_cast<size_t>(Sta[i]->state_nr * M), static_cast<size_t>(M));
				if (!(Seg[j]->freedom == "frozen" || chi == 0))
				{
					for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
					for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
					for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
				}
			}
		}
	}



	for (int j = 0; j < n_states; j++)
	{
		auto phi = std::span<const Real>(Seg[Sta[j]->mon_nr]->phi_state).subspan(static_cast<size_t>(Sta[j]->state_nr * M), static_cast<size_t>(M));
		for (int k = 0; k < n_seg; k++)
			if (Seg[k]->ns < 2)
			{
					chi = Sta[j]->chi[k];
					if (Seg[k]->freedom != "frozen") chi = chi / 2;


				auto phi_side = std::span<const Real>(Seg[k]->phi_side);
				if (chi != 0)
				{
					for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
					for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
					for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
				}
			}
		for (int l = 0; l < n_states; l++)
		{
				chi = Sta[j]->chi[n_seg + l] / 2;

			auto phi_side = std::span<const Real>(Seg[Sta[l]->mon_nr]->phi_side).subspan(static_cast<size_t>(Sta[l]->state_nr * M), static_cast<size_t>(M));
			if (chi != 0)
			{
				for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
				for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
				for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
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
		auto phi = std::span<const Real>(Mol[i]->phitot);
		std::copy(phi.begin(), phi.end(), temp.begin());
		for (int __i = 0; __i < M; ++__i) temp[__i] *= constant;
		for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
	}
	for (int __i = 0; __i < M; ++__i) F[__i] *= ksam[__i]; //clean up contributions in frozen sites.

std::fill(temp.begin(), temp.end(), 0.0);
	if (charged) {
		for (int __i = 0; __i < M; ++__i) temp[__i] += q[__i] * psi[__i];
		for (int __i = 0; __i < M; ++__i) temp[__i] *= 0.5;
		for (int __i = 0; __i < M; ++__i) F[__i] += temp[__i];
	}
	return FreeEnergy + lat->WeightedSum(F.data());
}

Real System::GetGrandPotential(void)
{ //Eqn 293
	NAMICS_DBG( "GetGrandPotential for system " << std::endl);
	int M = lat->M;
	auto& GP = GrandPotentialDensity;
	auto& temp = TEMP;
	auto& alpha = this->alpha;
	auto& EE = this->EE;
	auto& eps = this->eps;
	auto& q = this->q;
	auto& psi = this->psi;
	auto& ksam = KSAM;
	int n_mol = In->MolList.size();
	std::fill(GP.begin(), GP.end(), 0.0);

	for (int i = 0; i < n_mol; i++)
	{
		auto phi = std::span<const Real>(Mol[i]->phitot);
		Real phibulk = Mol[i]->phibulk;
		int N = Mol[i]->chainlength;
		std::copy(phi.begin(), phi.end(), temp.begin());
		for (int __i = 0; __i < M; ++__i) temp[__i] += -phibulk;
		for (int __i = 0; __i < M; ++__i) temp[__i] *= 1.0 / N; //GP has wrong sign. will be corrected at end of this routine;
		for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
	}

	for (int __i = 0; __i < M; ++__i) GP[__i] += alpha[__i];
	Real phibulkA;
	Real phibulkB;
	Real chi;
	int n_seg = In->MonList.size();
	int n_states = In->StateList.size();

	int n_mon = In->MonList.size();
	for (int i = 0; i < n_mon; i++) //if this is not done, in 3 gradients we have wrong results...
	{
		if (Seg[i]->ns < 2)
			lat->remove_bounds(Seg[i]->phi_side.data());
		else
		{
			for (int j = 0; j < Seg[i]->ns; j++)
			{
				lat->remove_bounds(Seg[i]->phi_side.data() + j * M);
			}
		}
	}


		for (int j = 0; j < n_seg; j++)
			if (Seg[j]->freedom != "frozen")
			{
				if (Seg[j]->ns < 2)
				{
				auto phi = std::span<const Real>(Seg[j]->phi);
				phibulkA = Seg[j]->phibulk;
					for (int k = 0; k < n_seg; k++)
					{
						if (Seg[k]->freedom != "frozen")
						{
						if (Seg[k]->ns < 2)
						{
								chi = Seg[j]->chi[k] / 2;

								auto phi_side = std::span<const Real>(Seg[k]->phi_side);
							phibulkB = Seg[k]->phibulk;
							if (chi != 0)
							{
								for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
								for (int __i = 0; __i < M; ++__i) temp[__i] += -phibulkA * phibulkB;
								for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
								for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
							}
						}
					}
				}

				for (int k = 0; k < n_states; k++)
				{
						chi = Seg[j]->chi[n_seg + k] / 2;

						auto phi_side = std::span<const Real>(Seg[Sta[k]->mon_nr]->phi_side).subspan(static_cast<size_t>(Sta[k]->state_nr * M), static_cast<size_t>(M));
					phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
					if (chi != 0)
					{
						for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
						for (int __i = 0; __i < M; ++__i) temp[__i] += -phibulkA * phibulkB;
						for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
						for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
					}
				}
			}
	}
	for (int j=0; j< n_seg; j++)
	{
		auto phi = std::span<const Real>(Seg[j]->phi);
		auto u_ext = std::span<const Real>(Seg[j]->u_ext);
		for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * u_ext[__i];
		for (int __i = 0; __i < M; ++__i) GP[__i] -= temp[__i];
	}


	for (int j = 0; j < n_states; j++)
	{
		auto phi = std::span<const Real>(Seg[Sta[j]->mon_nr]->phi_state).subspan(static_cast<size_t>(Sta[j]->state_nr * M), static_cast<size_t>(M));
		phibulkA = Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
			for (int k = 0; k < n_seg; k++)
				if (Seg[k]->freedom != "frozen")
			{
				if (Seg[k]->ns < 2)
				{
						chi = Sta[j]->chi[k] / 2;

						phibulkB = Seg[k]->phibulk;
					auto phi_side = std::span<const Real>(Seg[k]->phi_side);
					if (chi != 0)
					{
						for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
						for (int __i = 0; __i < M; ++__i) temp[__i] += -phibulkA * phibulkB;
						for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
						for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
					}
				}
			}
		for (int k = 0; k < n_states; k++)
		{
				chi = Sta[j]->chi[n_seg + k] / 2;

			auto phi_side = std::span<const Real>(Seg[Sta[k]->mon_nr]->phi_side).subspan(static_cast<size_t>(Sta[k]->state_nr * M), static_cast<size_t>(M));
			phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
			if (chi != 0)
			{
				for (int __i = 0; __i < M; ++__i) temp[__i] = phi[__i] * phi_side[__i];
				for (int __i = 0; __i < M; ++__i) temp[__i] += -phibulkA * phibulkB;
				for (int __i = 0; __i < M; ++__i) temp[__i] *= chi;
				for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
			}
		}
	}

	

	for (int __i = 0; __i < M; ++__i) GP[__i] *= -1.0; //correct the sign.

	std::fill(temp.begin(), temp.end(), 0.0);
if (charged) {
	for (int __i = 0; __i < M; ++__i) temp[__i] = EE[__i] * eps[__i];
	for (int __i = 0; __i < M; ++__i) temp[__i] *= -2.0;

	for (int __i = 0; __i < M; ++__i) temp[__i] += q[__i] * psi[__i];
	for (int __i = 0; __i < M; ++__i) temp[__i] *= -1.0/2.0;


		for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
		for (int __i = 0; __i < M; ++__i) GP[__i] *= ksam[__i];
		for (int __i = 0; __i < M; ++__i) temp[__i] = q[__i] * ksam[__i];
		for (int __i = 0; __i < M; ++__i) temp[__i] = q[__i] - temp[__i];
		for (int __i = 0; __i < M; ++__i) temp[__i] *= psi[__i];
		for (int __i = 0; __i < M; ++__i) temp[__i] *= 0.5;
	for (int __i = 0; __i < M; ++__i) GP[__i] += temp[__i];
}

	if (!charged) for (int __i = 0; __i < M; ++__i) GP[__i] *= ksam[__i]; //exclude solid and frozen sites from GP.

	return  lat->WeightedSum(GP.data());

}

bool System::CreateMu(int pos)
{
	NAMICS_DBG( "CreateMu for system " << std::endl);
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
		if (pos==M) Mu = std::log(NA * n / GN) + 1; else {
			Mu=std::log(Mol[i]->phitot[pos]) +1;
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
