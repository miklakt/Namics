#include "tools_host.h"
#include "input.h"
#include "lattice.h"
#include "LGrad1.h"
#include "LGrad2.h"
#include "LGrad3.h"
#include "LG1Planar.h"
#include "LG2Planar.h"
#include "molecule.h"
#include "mol_branched.h"
#include "mol_linear.h"
#include "namics.h"
#include "output.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "system.h"
#include "sfnewton.h"
#include "solve_scf.h"
#include <memory>

Real e = 1.60217e-19;
Real T = 298.15;
Real k_B = 1.38065e-23;
Real k_BT = k_B * T;
Real eps0 = 8.85418e-12;
Real PIE = 3.14159265;
//Used for command line switches
bool debug = false;

//Output when the user malforms input. Update when adding new command line switches.
void improperInput()
{
	std::cerr << "Improper usage: namics [-options] [filename]." << std::endl
			 << "Options available:" << std::endl;
	std::cerr << "-d Enables debugging mode." << std::endl;
}

int main(int argc, char *argv[])
{
	std::vector<std::string> args(argv, argv + argc);
	//Output error if no filename has been specified.
	if (argc == 1)
	{
		improperInput();
		return 1;
	}
	//Output error if user starts with a commandline switch. (Also catches combination with forgotten filename)
	if ((args.back())[0] == '-')
	{
		improperInput();
		return 1;
	}

	// If the specified filename has no extension: add the extension specified below.
	std::string extension = "in";
	std::ostringstream filename;
	filename << args.back();
	bool hasNoExtension = (filename.str().substr(filename.str().find_last_of(".") + 1) != extension);
	if (hasNoExtension)
		filename << "." << extension;

	//If the switch -d is given, enable debug. Add new switches by copying and replacing -d and debug = true.
	if (std::find(args.begin(), args.end(), "-d") != args.end())
	{
		debug = true;
	}

	int start = 0;
	int n_starts = 0;

	std::string METHOD = "";
	std::vector<Real> X;
	int MX = 0, MY = 0, MZ = 0;
	int fjc_old = 0;
	bool CHARGED = false;
	std::vector<std::string> MONLIST;
	std::vector<std::string> STATELIST;

	// Single ownership
	std::unique_ptr<Input> In;              // Inputs read from file
	std::unique_ptr<Lattice> Lat;
	std::unique_ptr<Lattice> lat_p;
	std::unique_ptr<Molecule> mol_p;
	std::unique_ptr<Solve_scf> New;         // Solver and iteration scheme
	std::unique_ptr<System> Sys;

	// Multi-instance collections
	std::unique_ptr<Output> Out;            // Output written to file
	std::vector<std::unique_ptr<Molecule>> Mol;
	std::vector<std::unique_ptr<Segment>> Seg;
	std::vector<std::unique_ptr<State>> Sta;
	std::vector<std::unique_ptr<Reaction>> Rea;

	In = std::make_unique<Input>(filename.str());
	if (In->Input_error)
	{
		return 0;
	}
	n_starts = In->GetNumStarts();
	if (n_starts == 0)
		n_starts++; // Default to 1 start..

	/******** This while loop basically contains the rest of main, initializes all classes and performs  ********/
	/******** calculations for a given number (start) of cycles ********/

	while (start < n_starts)
	{

		start++;
		if (!In->MakeLists(start)) return 0;
		std::cout << "Problem nr " << start << " out of " << n_starts << std::endl;

		/******** Class creation starts here ********/

		lat_p = std::make_unique<LGrad1>(*In, In->LatList[0]);
		if (!lat_p->CheckInput(start,true)) //-1 means that checkinput will stop when gradients and geometry are known.
		{
			return 0;
		} else
		{ int gradients=lat_p->gradients;
		  std::string geometry = lat_p->geometry;
		  bool success;
			lat_p.reset();
			switch (gradients) {
				case 1:
					if (geometry=="planar") {
						Lat = std::make_unique<LG1Planar>(*In,In->LatList[0]);
					} else {
						Lat = std::make_unique<LGrad1>(*In,In->LatList[0]);
					}

					break;
				case 2:
					if (geometry=="planar") {
						Lat = std::make_unique<LG2Planar>(*In,In->LatList[0]);
					} else {
						Lat = std::make_unique<LGrad2>(*In,In->LatList[0]);
					}
					break;
				case 3:
					Lat = std::make_unique<LGrad3>(*In,In->LatList[0]);
					break;
			}
			success=Lat->CheckInput(start,false);

			if (!success) return 0;
		}

		int n_seg = In->MonList.size();
		Seg.clear();
		Seg.reserve(n_seg);
		for (int i = 0; i < n_seg; i++) {
			Seg.push_back(std::make_unique<Segment>(In.get(), Lat.get(), In->MonList[i], i, n_seg));
		}
		//Create state class instance and check inputs
		int n_stat = In->StateList.size();
		Sta.clear();
		Sta.reserve(n_stat);
		for (int i = 0; i < n_stat; i++)
		{
			Sta.push_back(std::make_unique<State>(In.get(), Seg, In->StateList[i]));
		}

		for (int i = 0; i < n_seg; i++)
		{
			for (int k = 0; k < n_seg; k++)
			{
				Seg[i]->PutChiKEY(Seg[k]->name);
			}
			for (int k = 0; k < n_stat; k++)
			{
				Seg[i]->PutChiKEY(Sta[k]->name);
			}
			if (!Seg[i]->CheckInput(start))
				return 0;
		}
		for (int i = 0; i < n_stat; i++)
		{
			for (int k = 0; k < n_seg; k++)
			{
				Sta[i]->PutChiKEY(Seg[k]->name);
			}
			for (int k = 0; k < n_stat; k++)
			{
				Sta[i]->PutChiKEY(Sta[k]->name);
			}
			if (!Sta[i]->CheckInput(start))
				return 0;
		}

		//Create reaction class instance and check inputs
		int n_rea = In->ReactionList.size();
		Rea.clear();
		Rea.reserve(n_rea);
		for (int i = 0; i < n_rea; i++)
		{
			Rea.push_back(std::make_unique<Reaction>(In.get(), Seg, Sta, In->ReactionList[i]));
			if (!Rea[i]->CheckInput(start))
				return 0;
		}

		int n_mol = In->MolList.size();
		Mol.clear();
		Mol.reserve(n_mol);
		for (int i = 0; i < n_mol; i++)
		{
			mol_p = std::make_unique<Molecule>(In.get(), Lat.get(), Seg, In->MolList[i]);
			if (!mol_p->CheckInput(start,true)) //'true' here means that checkinput can stop wehn Moltype and freedom are known.
			{
				return 0;
			} else {
				if (mol_p->MolType == monomer) {
					Mol.push_back(std::make_unique<Molecule>(In.get(), Lat.get(), Seg, In->MolList[i]));
				} else if (mol_p->MolType == linear) {
					Mol.push_back(std::make_unique<mol_linear>(In.get(), Lat.get(), Seg, In->MolList[i]));
				} else {
					Mol.push_back(std::make_unique<mol_branched>(In.get(), Lat.get(), Seg, In->MolList[i]));
				}
				mol_p.reset();
				if (!Mol[i]->CheckInput(start,false)) return 0;
			}
		}

		Sys = std::make_unique<System>(In.get(), Lat.get(), Seg, Sta, Rea, Mol, In->SysList[0]);
		if (!Sys->CheckInput(start)) return 0;
		if (!Sys->CheckChi_values(n_seg))return 0;


		New = std::make_unique<Solve_scf>(In.get(), Lat.get(), Seg, Sta, Rea, Mol, Sys.get(), In->NewtonList[0]);
		if (!New->CheckInput(start)) return 0;



		//Guesses geometry
		if (Sys->initial_guess == "file")
		{
			MONLIST.clear();
			STATELIST.clear();
			if (!io::ReadInitialGuess(Sys->guess_inputfile, std::span<Real>{}, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old, 0))
			{
				return 1;
			}
			int nummon = MONLIST.size();
			int numstate = STATELIST.size();
			int m;
			if (MY == 0) m = MX + 2*fjc_old;
			else {
				if (MZ == 0) {
					m = (MX + 2*fjc_old) * (MY + 2*fjc_old);
				} else {
					m = (MX + 2*fjc_old) * (MY + 2*fjc_old) * (MZ + 2*fjc_old);
				}
			}
			int IV = (nummon + numstate) * m;

			if (CHARGED)
				IV += m;
			X.resize(IV);
			MONLIST.clear();
			STATELIST.clear();
			if (!io::ReadInitialGuess(Sys->guess_inputfile, std::span<Real>(X), METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old, 1)) {
				return 1;
			}
		}

		int IV_new=0;
		int substart = 0;
		int subloop = 0;
		int mon_length;
		int state_length;
		// Prepare and create output class instance.
		Out.reset();
		if (In->OutputList.empty()) {
			std::cout << "Warning: no output defined!" << std::endl;
		} else {
			Out = std::make_unique<Output>(In.get(), Lat.get(), Seg, Sta, Rea, Mol, Sys.get(), New.get(), "json");
			if (!Out->CheckInput(start)) {
				std::cout << "input_error in output " << std::endl;
				return 0;
			}
		}

		while (subloop <= substart)
		{
			Sys->MakeItsLists();

				New->AllocateMemory();

				if (Sys->initial_guess != "none" && !X.empty())
				New->Guess(X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);

				if (!New->Solve(true)) return 1;

				if (Sys->initial_guess == "previous_result" || Sys->initial_guess == "file") {
					if (New->iv == IV_new) {
						std::copy_n(New->xx.begin(), IV_new, X.begin());
					} else {
						IV_new=New->iv;
						X.resize(IV_new);
						std::copy_n(New->xx.begin(), IV_new, X.begin());
						MX=Lat->MX;
						MY=Lat->MY;
						MZ=Lat->MZ;
						fjc_old=Lat->fjc;
						mon_length = Sys->ItMonList.size();
						state_length = Sys->ItStateList.size();
						MONLIST.clear();
						STATELIST.clear();
						for (int i = 0; i < mon_length; i++)
						{
							MONLIST.push_back(Seg[Sys->ItMonList[i]]->name);
						}
						for (int i = 0; i < state_length; i++)
						{
							STATELIST.push_back(Sta[Sys->ItStateList[i]]->name);
						}
					}
				}
				New->PushOutput();

				if (Out) Out->WriteOutput(subloop);


				subloop++;
			}

		if (Sys->initial_guess == "previous_result"|| Sys->initial_guess == "file")
		{
			METHOD = New->SCF_method;
			MX = Lat->MX;
			MY = Lat->MY;
			MZ = Lat->MZ;
			CHARGED = Sys->charged;
			IV_new = New->iv;
			X.resize(IV_new);
			std::copy_n(New->xx.begin(), IV_new, X.begin());
			fjc_old = Lat->fjc;
			mon_length = Sys->ItMonList.size();
			state_length = Sys->ItStateList.size();
			MONLIST.clear();
			STATELIST.clear();
			for (int i = 0; i < mon_length; i++)
			{
				MONLIST.push_back(Seg[Sys->ItMonList[i]]->name);
			}
			for (int i = 0; i < state_length; i++)
			{
				STATELIST.push_back(Sta[Sys->ItStateList[i]]->name);
			}
		}
		/******** Clear all class instances ********/

		lat_p.reset();
		mol_p.reset();

		Out.reset();
		New.reset();
		Sys.reset();
		Mol.clear();
		Seg.clear();
		Sta.clear();
		Rea.clear();
		Lat.reset();
	} //loop over starts.
	In.reset();
	return 0;
}
