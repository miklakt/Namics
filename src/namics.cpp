#include "tools.h"
#include "input.h"
#include "input_preprocessor.h"
#include "io_utils.h"
#include "lattice.h"
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
#include <filesystem>
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

	std::filesystem::path filename(args.back());
	if (!filename.has_extension()) filename += ".in";

	//If the switch -d is given, enable debug. Add new switches by copying and replacing -d and debug = true.
	if (std::find(args.begin(), args.end(), "-d") != args.end())
	{
		debug = true;
	}

	int start = 0;
	int n_starts = 0;

	std::vector<Real> X;
	int MX = 0, MY = 0, MZ = 0;
	int fjc_old = 0;
	bool CHARGED = false;
	std::vector<std::string> MONLIST;
	std::vector<std::string> STATELIST;

	// Single ownership
	std::unique_ptr<Input> In;              // Inputs read from file
	std::unique_ptr<Lattice> Lat;
	std::unique_ptr<Molecule> mol_p;
	std::unique_ptr<Solve_scf> New;         // Solver and iteration scheme
	std::unique_ptr<System> Sys;

	// Multi-instance collections
	std::unique_ptr<Output> Out;            // Output written to file
	std::vector<std::unique_ptr<Molecule>> Mol;
	std::vector<std::unique_ptr<Segment>> Seg;
	std::vector<std::unique_ptr<State>> Sta;
	std::vector<std::unique_ptr<Reaction>> Rea;

	std::string json_input_path;
	if (!io::input::PrepareInputFile(filename.string(), json_input_path)) {
		return 0;
	}

	In = std::make_unique<Input>(json_input_path);
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

		Lat = lattice_factory::CreateChecked(*In, In->LatList[0], start);
		if (!Lat) return 0;

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
			if (!Seg[i]->CheckInput(start))
				return 0;
		}
		for (int i = 0; i < n_stat; i++)
		{
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
			for (const int mon_index : Sys->ItMonList) MONLIST.push_back(Seg[mon_index]->name);
			for (const int state_index : Sys->ItStateList) STATELIST.push_back(Sta[state_index]->name);
			CHARGED = Sys->charged;
			const int IV = static_cast<int>(MONLIST.size() + STATELIST.size() + (CHARGED ? 1 : 0)) * Lat->M;
			X.resize(IV);
			Lat->AllocateMemory();
			auto pad_profile = [&](std::vector<Real>& values) {
				const int trimmed_size = Lat->MX * (Lat->gradients >= 2 ? Lat->MY : 1) * (Lat->gradients >= 3 ? Lat->MZ : 1);
				if (static_cast<int>(values.size()) != trimmed_size) return false;
				std::vector<Real> padded(static_cast<size_t>(Lat->M), 0);
				size_t offset = 0;
				if (Lat->gradients == 1) {
					for (int x = Lat->fjc; x < Lat->MX + Lat->fjc; ++x) padded[x] = values[offset++];
				} else if (Lat->gradients == 2) {
					for (int x = Lat->fjc; x < Lat->MX + Lat->fjc; ++x) {
						for (int y = Lat->fjc; y < Lat->MY + Lat->fjc; ++y) padded[Lat->P(x, y)] = values[offset++];
					}
				} else {
					for (int x = Lat->fjc; x < Lat->MX + Lat->fjc; ++x) {
						for (int y = Lat->fjc; y < Lat->MY + Lat->fjc; ++y) {
							for (int z = Lat->fjc; z < Lat->MZ + Lat->fjc; ++z) padded[Lat->P(x, y, z)] = values[offset++];
						}
					}
				}
				Lat->set_bounds(padded.data());
				values.swap(padded);
				return true;
			};
			auto resolve_profile = [&](const io::detail::json& document, const std::string& key, std::vector<Real>& values) {
				const io::detail::json* source = io::detail::FindLastProblemObject(document);
				if (source == nullptr && document.is_object()) source = &document;
				if (source == nullptr) return false;
				if (key == "psi") {
					const auto sys_section = source->find("sys");
					if (sys_section == source->end() || !sys_section->is_object()) return false;
					const io::detail::json* profile = io::detail::FindProfileInOutputSection(*sys_section, Sys->name, "psi");
					if (profile == nullptr) {
						for (auto sys = sys_section->begin(); sys != sys_section->end() && profile == nullptr; ++sys) {
							profile = io::detail::FindProfileInOutputSection(*sys_section, sys.key(), "psi");
						}
					}
					return profile != nullptr && io::detail::ReadOutputProfileArray(*profile, Lat->M, values, pad_profile);
				}
				const bool is_state = key.rfind("state:", 0) == 0;
				const std::string name = key.substr(is_state ? 6 : 4);
				const auto& section = is_state ? source->find("state") : source->find("mon");
				if (section == source->end() || !section->is_object()) return false;
				const io::detail::json* profile = io::detail::FindProfileInOutputSection(*section, name, "u");
				if (profile == nullptr) {
					if (is_state) {
						int state_index = -1;
						for (size_t i = 0; i < Sta.size(); ++i) if (Sta[i]->name == name) { state_index = static_cast<int>(i); break; }
						if (state_index < 0) return false;
						for (size_t i = 0; i < Sta.size() && profile == nullptr; ++i) if (Sta[i]->chi == Sta[state_index]->chi) profile = io::detail::FindProfileInOutputSection(*section, Sta[i]->name, "u");
					} else {
						int mon_index = -1;
						for (size_t i = 0; i < Seg.size(); ++i) if (Seg[i]->name == name) { mon_index = static_cast<int>(i); break; }
						if (mon_index < 0) return false;
						for (size_t i = 0; i < Seg.size() && profile == nullptr; ++i) if (Seg[i]->chi == Seg[mon_index]->chi) profile = io::detail::FindProfileInOutputSection(*section, Seg[i]->name, "u");
					}
				}
				return profile != nullptr && io::detail::ReadOutputProfileArray(*profile, Lat->M, values, pad_profile);
			};
			if (!io::ReadInitialGuess(Sys->guess_inputfile, std::span<Real>(X), MONLIST, STATELIST, CHARGED, resolve_profile)) {
				return 1;
			}
		}

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

		Sys->MakeItsLists();
		New->AllocateMemory();

		if (!X.empty()) {
			if (Sys->initial_guess == "file") {
				if (static_cast<int>(X.size()) != New->iv) {
					std::cout << "Input initial_guess size does not match the current system." << std::endl;
					return 1;
				}
				std::copy_n(X.begin(), New->iv, New->xx.begin());
			} else if (Sys->initial_guess != "none") {
				New->Guess(X, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
			}
		}

		if (!New->Solve(true)) return 1;
		New->PushOutput();

		if (Out) Out->WriteOutput(0);

		if (Sys->initial_guess == "previous_result")
		{
			MX = Lat->MX;
			MY = Lat->MY;
			MZ = Lat->MZ;
			CHARGED = Sys->charged;
			X.resize(New->iv);
			std::copy_n(New->xx.begin(), New->iv, X.begin());
			fjc_old = Lat->fjc;
			MONLIST.clear();
			STATELIST.clear();
			for (const int mon_index : Sys->ItMonList) MONLIST.push_back(Seg[mon_index]->name);
			for (const int state_index : Sys->ItStateList) STATELIST.push_back(Sta[state_index]->name);
		}
		/******** Clear all class instances ********/

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
