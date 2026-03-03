#define MAINxH
#include "tools_host.h"
#include "input.h"
#include "lattice.h"
//#include "lat_preview.h"
//#include "mol_preview.h"
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

string version = "2.2.2.2.2.1.1";
// meaning:
// newton version number =2
// system version number =2
// lattice version number =2
// molecule version number =2
// segment version number =2
// alias version number =1
// output version number =1
Real check = 0.4534345;
Real e = 1.60217e-19;
Real T = 298.15;
Real k_B = 1.38065e-23;
Real k_BT = k_B * T;
Real eps0 = 8.85418e-12;
Real PIE = 3.14159265;
int DEBUG_BREAK = 1;
//Used for command line switches
bool debug = false;

//Output when the user malforms input. Update when adding new command line switches.
void improperInput()
{
	cerr << "Improper usage: namics [-options] [filename]." << endl
			 << "Options available:" << endl;
	cerr << "-d Enables debugging mode." << endl;
}

int main(int argc, char *argv[])
{
	vector<string> args(argv, argv + argc);
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
	string extension = "in";
	ostringstream filename;
	filename << args.back();
	bool hasNoExtension = (filename.str().substr(filename.str().find_last_of(".") + 1) != extension);
	if (hasNoExtension)
		filename << "." << extension;

	//If the switch -d is given, enable debug. Add new switches by copying and replacing -d and debug = true.
	if (find(args.begin(), args.end(), "-d") != args.end())
	{
		debug = true;
	}

	int start = 0;
	int n_starts = 0;

	string final_guess;
	string METHOD = "";
	Real *X = NULL;
	int MX = 0, MY = 0, MZ = 0;
	int fjc_old = 0;
	bool CHARGED = false;
	vector<string> MONLIST;
	vector<string> STATELIST;

	// Single ownership
	unique_ptr<Input> In;              // Inputs read from file
	unique_ptr<Lattice> Lat;
	unique_ptr<Lattice> lat_p;
	unique_ptr<Molecule> mol_p;
	unique_ptr<Solve_scf> New;         // Solver and iteration scheme
	unique_ptr<System> Sys;

	// Multi-instance collections
	vector<Output *> Out;              // Outputs written to file
	vector<Molecule *> Mol;            // Properties of entire molecule
	vector<Segment *> Seg;             // Properties of molecule segments
	vector<State *> Sta;
	vector<Reaction *> Rea;

	// Create input class instance and handle errors(reference above)
	In = make_unique<Input>(filename.str());
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
		cout << "Problem nr " << start << " out of " << n_starts << endl;

		/******** Class creation starts here ********/

		// Create lattice class instance and check inputs (reference above)
		lat_p = make_unique<LGrad1>(*In, In->LatList[0]);

		//Lat->outputtest();
		if (!lat_p->CheckInput(start,true)) //-1 means that checkinput will stop when gradients and geometry are known.
		{
			return 0;
		} else
		{ int gradients=lat_p->gradients;
		  string geometry = lat_p->geometry;
		  bool success;
			lat_p.reset();
			switch (gradients) {
				case 1:
					if (geometry=="planar") {
						Lat = make_unique<LG1Planar>(*In,In->LatList[0]);
					} else {
						Lat = make_unique<LGrad1>(*In,In->LatList[0]);
					}

					break;
				case 2:
					if (geometry=="planar") {
						Lat = make_unique<LG2Planar>(*In,In->LatList[0]);
					} else {
						Lat = make_unique<LGrad2>(*In,In->LatList[0]);
					}
					break;
				case 3:
					Lat = make_unique<LGrad3>(*In,In->LatList[0]);
					break;
				default :
					break;

			}
			success=Lat->CheckInput(start,false);

			if (!success) return 0;
		}

		// Create segment class instance and check inputs (reference above)
		int n_seg = In->MonList.size();
		for (int i = 0; i < n_seg; i++) {
			Seg.push_back(new Segment(In.get(), Lat.get(), In->MonList[i], i, n_seg));
		}
		//Create state class instance and check inputs
		int n_stat = In->StateList.size();
		for (int i = 0; i < n_stat; i++)
			Sta.push_back(new State(In.get(), Seg, In->StateList[i]));

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
		for (int i = 0; i < n_rea; i++)
		{
			Rea.push_back(new Reaction(In.get(), Seg, Sta, In->ReactionList[i]));
			if (!Rea[i]->CheckInput(start))
				return 0;
		}

		// Create segment class instance and check inputs (reference above)
		int n_mol = In->MolList.size();
		for (int i = 0; i < n_mol; i++)
		{
			mol_p = make_unique<Molecule>(In.get(), Lat.get(), Seg, In->MolList[i]);
			if (!mol_p->CheckInput(start,true)) //'true' here means that checkinput can stop wehn Moltype and freedom are known.
			{
				return 0;
			} else {
				switch (mol_p->MolType) {

				 	case monomer:
						Mol.push_back(new Molecule(In.get(), Lat.get(), Seg, In->MolList[i]));
						break;
						case water:
							Mol.push_back(new Molecule(In.get(), Lat.get(), Seg, In->MolList[i]));
							break;
					case linear:
						if (mol_p->freedom=="clamped") {
							cout << "Unsupported molecule freedom 'clamped' for mol '" << In->MolList[i] << "'." << endl;
							return 0;
						}
						Mol.push_back(new mol_linear(In.get(), Lat.get(), Seg, In->MolList[i]));
						break;
					case branched:
						Mol.push_back(new mol_branched(In.get(), Lat.get(), Seg, In->MolList[i]));
						break;
					case dendrimer:
					case asym_dendrimer:
					case comb:
						cout << "Unsupported molecule architecture for mol '" << In->MolList[i]
						     << "'. Only monomer, linear, branched and water are supported." << endl;
						return 0;
						break;
					default:
						cout <<"Unknown MolType " << endl;
						break;

				}
				mol_p.reset();
				if (!Mol[i]->CheckInput(start,false)) return 0;
			}
		}

		// Create system class instance and check inputs (reference above)
		Sys = make_unique<System>(In.get(), Lat.get(), Seg, Sta, Rea, Mol, In->SysList[0]);
		if (!Sys->CheckInput(start)) return 0;
		if (!Sys->CheckChi_values(n_seg))return 0;


		// Create newton class instance and check inputs (reference above)
		New = make_unique<Solve_scf>(In.get(), Lat.get(), Seg, Sta, Rea, Mol, Sys.get(), In->NewtonList[0]);
		if (!New->CheckInput(start)) return 0;



		//Guesses geometry
		if (Sys->initial_guess == "file")
		{
			MONLIST.clear();
			STATELIST.clear();
			if (!Lat->ReadGuess(Sys->guess_inputfile, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old, 0))
			{
				// last argument 0 is to first checkout sizes of system.
				return 0;
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
			if (start > 0) {
				free(X);
				X = (Real *)malloc(IV * sizeof(Real));
			}
			MONLIST.clear();
			STATELIST.clear();
			Lat->ReadGuess(Sys->guess_inputfile, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old, 1);
			// last argument 1 is to read guess in X.
		}

		int IV_new=0;
		int substart = 0;
		int subloop = 0;
		int n_out = 0;
		int mon_length;
		int state_length;
			// Prepare, catch errors for output class creation
			n_out = In->OutputList.size();
			if (n_out == 0)
				cout << "Warning: no output defined!" << endl;

			// Create output class instance and check inputs (reference above)
			for (int ii = 0; ii < n_out; ii++)
			{
				Out.push_back(new Output(In.get(), Lat.get(), Seg, Sta, Rea, Mol, Sys.get(), New.get(), In->OutputList[ii], ii, n_out));
				if (!Out[ii]->CheckInput(start))
				{
					cout << "input_error in output " << endl;
					return 0;
				}
			}

			while (subloop <= substart)
			{
				Sys->MakeItsLists();

				New->AllocateMemory();
				//} else

				if (Sys->initial_guess != "none")
				New->Guess(X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);

				New->Solve(true);

				if (Sys->initial_guess == "previous_result" || Sys->initial_guess == "file") {
					if (New->iv == IV_new) {
						std::copy_n(New->xx, IV_new, X);
					} else {
						if (X!=NULL) free(X);
						IV_new=New->iv;
						X = (Real *)malloc(IV_new * sizeof(Real));
						std::copy_n(New->xx, IV_new, X);
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

				for (int ii = 0; ii < n_out; ii++)
				{
					Out[ii]->WriteOutput(subloop);
				}
				if (Sys->final_guess == "file")
				{
					MONLIST.clear();
					STATELIST.clear();
					int mon_length = Sys->ItMonList.size();
					int state_length = Sys->ItStateList.size();
					for (int i = 0; i < mon_length; i++)
					{
						MONLIST.push_back(Seg[Sys->ItMonList[i]]->name);
					}
					for (int i = 0; i < state_length; i++)
					{
						STATELIST.push_back(Sta[Sys->ItStateList[i]]->name);
					}
					Lat->StoreGuess(Sys->guess_outputfile, New->xx, New->SCF_method, MONLIST, STATELIST, Sys->charged, start);
				}


				subloop++;
			}

		for (auto all_segments : Seg)
			all_segments->prepared = false;

		if (Sys->initial_guess == "previous_result"|| Sys->initial_guess == "file")
		{
			METHOD = New->SCF_method; //check this..
			MX = Lat->MX;
			MY = Lat->MY;
			MZ = Lat->MZ;
			CHARGED = Sys->charged;
			IV_new = New->iv; //check this
			if (start > 1 || (start == 1 && Sys->initial_guess == "file"))
				free(X);
			X = (Real *)malloc(IV_new * sizeof(Real));
			std::copy_n(New->xx, IV_new, X);
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
		if (Sys->final_guess == "file")
		{
			MONLIST.clear();
			STATELIST.clear();
			int mon_length = Sys->ItMonList.size();
			int state_length = Sys->ItStateList.size();
			for (int i = 0; i < mon_length; i++)
			{
				MONLIST.push_back(Seg[Sys->ItMonList[i]]->name);
			}
			for (int i = 0; i < state_length; i++)
			{
				STATELIST.push_back(Sta[Sys->ItStateList[i]]->name);
			}
			Lat->StoreGuess(Sys->guess_outputfile, New->xx, New->SCF_method, MONLIST, STATELIST, Sys->charged, start);
		}

		/******** Clear all class instances ********/

		lat_p.reset();
		mol_p.reset();

		for (int i = 0; i < n_out; i++)
			delete Out[i];
		Out.clear();
		New.reset();
		Sys.reset();
		for (int i = 0; i < n_mol; i++)
			delete Mol[i];
		Mol.clear();
		for (int i = 0; i < n_seg; i++)
			delete Seg[i];
		Seg.clear();
		for (int i = 0; i < n_stat; i++)
			delete Sta[i];
		Sta.clear();
		for (int i = 0; i < n_rea; i++)
			delete Rea[i];
		Rea.clear();
		Lat.reset();
	} //loop over starts.
	free(X);
	In.reset();
	return 0;
}
