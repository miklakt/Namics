#include "lattice.h"
#include "LGrad1.h"
#include "LGrad2.h"
#include "LGrad3.h"

#include <string_view>

namespace {

template <typename VecFn1, typename VecFn2, typename VecFn3>
void AssignBoundary(std::string_view bc,
                    int fjc,
                    int& single,
                    std::vector<int>& vector,
                    int mirror_single,
                    int periodic_single,
                    int surface_single,
                    VecFn1 mirror_vec,
                    VecFn2 periodic_vec,
                    VecFn3 surface_vec) {
	if (fjc == 1) {
		if (bc == "mirror") single = mirror_single;
		else if (bc == "periodic") single = periodic_single;
		else if (bc == "surface") single = surface_single;
		return;
	}
	if (bc == "mirror") {
		for (int k = 0; k < fjc; ++k) vector[k] = mirror_vec(k);
	} else if (bc == "periodic") {
		for (int k = 0; k < fjc; ++k) vector[k] = periodic_vec(k);
	} else if (bc == "surface") {
		for (int k = 0; k < fjc; ++k) vector[k] = surface_vec(k);
	}
}

template <class T>
std::unique_ptr<Lattice> TryCreate(const Input& in, const std::string& name, int start, const LatticeSelection& selection) {
	if (!T::Matches(selection)) return nullptr;
	auto lattice = std::make_unique<T>(in, name);
	if (!lattice->CheckInput(start)) return nullptr;
	return lattice;
}

} // namespace

Lattice::Lattice(const Input& In_,const std::string& name_) :
	BC(6) // boundary condition slots: lower/upper for x, y, z
{ //this file contains switch (gradients). In this way we keep all the lattice issues in one file!
NAMICS_DBG("Lattice constructor" << std::endl);	In=&In_; name=name_;
	sub_box_on = 0;
	all_lattice = false;
	ignore_sites=false;
	fcc_sites=false;
	stencil_full=false;
	gradients=1;
	fjc=1;
	MX=MY=MZ=0;
	offset_first_layer=0;
	subl=0;
}

namespace lattice_factory {

std::unique_ptr<Lattice> CreateChecked(const Input& in, const std::string& name, int start) {
	LatticeSelection selection;
	const auto& parameters = in.Parameters("lat", name, start);
	const std::vector<std::string> keys = {
		"gradients", "n_layers", "offset_first_layer", "geometry",
		"n_layers_x", "n_layers_y", "n_layers_z",
		"lowerbound", "upperbound",
		"lowerbound_x", "upperbound_x",
		"lowerbound_y", "upperbound_y",
		"lowerbound_z", "upperbound_z",
		"bondlength", "ignore_site_fraction", "fcc_site_fraction",
		"lattice_type", "stencil_full", "FJC_choices", "b/l"
	};
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		std::cout << "lat property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : keys) std::cout << item << std::endl;
		return nullptr;
	}
	try {
		selection.gradients = parameters.value("gradients", 1);
		if (selection.gradients < 1 || selection.gradients > 3) {
			std::cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << std::endl;
			selection.gradients = 1;
		}

		selection.geometry = parameters.value("geometry", std::string{"planar"});
		if (selection.geometry != "spherical" && selection.geometry != "cylindrical" && selection.geometry != "flat" && selection.geometry != "planar") {
			std::cout << "In lattice input for 'geometry' not recognized." << std::endl;
			return nullptr;
		}
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in lat '" << name << "': " << error.what() << std::endl;
		return nullptr;
	}
	if (selection.geometry == "flat") selection.geometry = "planar";

	if (auto lattice = TryCreate<LGrad1>(in, name, start, selection)) return lattice;
	if (auto lattice = TryCreate<LGrad2>(in, name, start, selection)) return lattice;
	if (auto lattice = TryCreate<LGrad3>(in, name, start, selection)) return lattice;

	std::cout << "No lattice implementation matches gradients=" << selection.gradients
	          << " and geometry=" << selection.geometry << std::endl;
	return nullptr;
}

} // namespace lattice_factory

void Lattice::DeAllocateMemory(void) {
NAMICS_DBG("DeAllocateMemory in lat " << std::endl);	if (!all_lattice) return;
	all_lattice=false;
	H.clear();
	B_X1.clear();
	B_Y1.clear();
	B_Z1.clear();
	B_XM.clear();
	B_YM.clear();
	B_ZM.clear();
	L.clear();
	lambda0.clear();
	fcc_lambda0.clear();
	lambda_1.clear();
	fcc_lambda_1.clear();
	lambda1.clear();
	fcc_lambda1.clear();
	LAMBDA.clear();
	X.clear();
}


void Lattice::AllocateMemory(void) {
NAMICS_DBG("AllocateMemory in lat " << std::endl);
	DeAllocateMemory();
	all_lattice=true;
	PutM();
	if (fjc>1) {
		B_X1.assign(fjc, 0);
		B_Y1.assign(fjc, 0);
		B_Z1.assign(fjc, 0);
		B_XM.assign(fjc, 0);
		B_YM.assign(fjc, 0);
		B_ZM.assign(fjc, 0);
	}

	switch (gradients) {
		case 3:
			AssignBoundary(BC[2], fjc, BZ1, B_Z1, 1, MZ, 0,
				[&](int k) { return 2 * fjc - 1 - k; },
				[&](int k) { return MZ + k; },
				[&](int k) { return k; });
			AssignBoundary(BC[5], fjc, BZM, B_ZM, MZ, 1, MZ + 1,
				[&](int k) { return MZ + fjc - k - 1; },
				[&](int k) { return fjc + k; },
				[&](int k) { return MZ + fjc + k; });

			//Fall through
		case 2:
			AssignBoundary(BC[1], fjc, BY1, B_Y1, 1, MY, 0,
				[&](int k) { return 2 * fjc - 1 - k; },
				[&](int k) { return MY + k; },
				[&](int k) { return k; });
			AssignBoundary(BC[4], fjc, BYM, B_YM, MY, 1, MY + 1,
				[&](int k) { return MY + fjc - k - 1; },
				[&](int k) { return fjc + k; },
				[&](int k) { return MY + fjc + k; });

			//Fall through
		case 1:
			AssignBoundary(BC[0], fjc, BX1, B_X1, 1, MX, 0,
				[&](int k) { return 2 * fjc - 1 - k; },
				[&](int k) { return MX + k; },
				[&](int k) { return k; });
			AssignBoundary(BC[3], fjc, BXM, B_XM, MX, 1, MX + 1,
				[&](int k) { return MX + fjc - k - 1; },
				[&](int k) { return fjc + k; },
				[&](int k) { return MX + fjc + k; });

	}
	if (fcc_sites) {
		fcc_lambda_1.assign(M, 0);
		fcc_lambda1.assign(M, 0);
		fcc_lambda0.assign(M, 0);
	}

	if (fjc==1) {
		if (gradients<3) {
		L.assign(M, 0);
		lambda_1.assign(M, 0);
		lambda1.assign(M, 0);
		lambda0.assign(M, 0);
		}
	} else {
		L.assign(M, 0);
		LAMBDA.assign(FJC * M, 0);
	}
	H.assign(M, 0);
	X.assign(M, 0);
	ComputeLambdas();
}



Lattice::~Lattice() {
NAMICS_DBG("lattice destructor " << std::endl); DeAllocateMemory();

}


int Lattice::P(int x, int y, int z) {
	return x*JX+y*JY+z;//(x+fjc-1)*JX + (y+fjc-1)*JY +(z+fjc-1);
}
int Lattice::P(int x, int y) {
	return x*JX+y; //(x+fjc-1)*JX +(y+fjc-1);
}

int Lattice::P(int x) {
	return x; //x+fjc-1;
}

bool Lattice::AssignChoice(const std::string& value, std::string& target, std::initializer_list<const char*> allowed, const char* error) const {
	for (const char* option : allowed) {
		if (value != option) continue;
		target = value;
		return true;
	}
	std::cout << error << std::endl;
	return false;
}

bool Lattice::ReadBoundaryCondition(const ParameterStore& parameters, const char* key, int slot, std::initializer_list<const char*> allowed, const char* error, const char* fallback) {
	const std::string value = parameters.value(key, std::string{});
	if (value.empty()) {
		BC[slot] = fallback;
		return true;
	}
	return AssignChoice(value, BC[slot], allowed, error);
}

bool Lattice::ReadScaledDimension(const ParameterStore& parameters, const char* key, int& value, int min_value, const char* missing_message, const char* bounds_message) {
	if (!parameters.contains(key)) {
		std::cout << missing_message << std::endl;
		return false;
	}
	value = parameters.at(key).get<int>();
	if (value < min_value || value > 1e6) {
		std::cout << bounds_message << std::endl;
		return false;
	}
	value *= fjc;
	return true;
}

void Lattice::ReadOffsetFirstLayer(const ParameterStore& parameters) {
	offset_first_layer = parameters.value("offset_first_layer", 0.0);
	if (offset_first_layer < 0) {
		std::cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << std::endl;
		offset_first_layer = 0;
	}
	offset_first_layer *= fjc;
}

bool Lattice::RejectParameters(const ParameterStore& parameters, std::initializer_list<std::pair<const char*, const char*>> rejected) const {
	bool success = true;
	for (const auto& [key, message] : rejected) {
		if (!parameters.contains(key)) continue;
		std::cout << message << std::endl;
		success = false;
	}
	return success;
}

bool Lattice::RejectAxisBoundsIn1D(const ParameterStore& parameters) const {
	return RejectParameters(parameters, {
		{"lowerbound_x", "lowerbound_x is not allowed in 1-gradient calculations"},
		{"lowerbound_y", "lowerbound_y is not allowed in 1-gradient calculations"},
		{"lowerbound_z", "lowerbound_z is not allowed in 1-gradient calculations"},
		{"upperbound_x", "upperbound_x is not allowed in 1-gradient calculations"},
		{"upperbound_y", "upperbound_y is not allowed in 1-gradient calculations"},
		{"upperbound_z", "upperbound_z is not allowed in 1-gradient calculations"},
	});
}

bool Lattice::RejectScalarBoundsInMultiD(const ParameterStore& parameters) const {
	return RejectParameters(parameters, {
		{"upperbound", "upperbound is only allowed in 1-gradient calculations"},
		{"lowerbound", "lowerbound is only allowed in 1-gradient calculations"},
	});
}

bool Lattice::RejectZBoundsIn2D(const ParameterStore& parameters) const {
	return RejectParameters(parameters, {
		{"lowerbound_z", "lowerbound_z is not allowed in 2-gradient calculations"},
		{"upperbound_z", "upperbound_z is not allowed in 2-gradient calculations"},
	});
}

bool Lattice::CheckPeriodicPair(int lower_slot, int upper_slot, const char* message) const {
	if (BC[lower_slot] != "periodic" && BC[upper_slot] != "periodic") return true;
	if (BC[lower_slot] == BC[upper_slot]) return true;
	std::cout << message << std::endl;
	return false;
}

bool Lattice::CheckInput(int start) {
NAMICS_DBG("CheckInput in lattice " << std::endl);	bool success=true;
	mx.push_back(0); my.push_back(0); mz.push_back(0); jx.push_back(0); jy.push_back(0); m.push_back(0); n_box.push_back(0);
	const auto& parameters = In->Parameters("lat", name, start);

	try {
		FJC=3;	fjc=1;
		if (parameters.contains("FJC_choices")) {
			FJC = parameters.at("FJC_choices").get<int>();
			if ((FJC-3) %2 != 0) {
				std::cout << "FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ...." <<std::endl;
				success=false;
			}
			fjc=(FJC-1)/2;
		}

		if (parameters.contains("b/l")) {
			int fjc_new = parameters.at("b/l").get<int>();
			if (fjc_new <1 ) {
				std::cout << "b/l should be a positive integer: 1, 2, 3, ...." <<std::endl;
				success=false;
			}
			if (parameters.contains("FJC_choices") && fjc_new !=fjc) {
				std::cout <<"You have set both 'FJC_choices' and 'b/l', but their values are not consistent with each other."<<std::endl;
				if (fjc_new<fjc && fjc_new >0) {
					std::cout <<"The value of 'b/l' is used, and that of FJC_choices is rejected." << std::endl;
				} else {
					if (fjc_new<1) success=true;
					std::cout <<"The value of 'FJC_choices' is used, and that of b/l is rejected." << std::endl;
					fjc_new=fjc;
				}
			}
			fjc=fjc_new;
			FJC=2*fjc+1;
		}

		bond_length=0;
		if (parameters.contains("bondlength")) {
			bond_length = parameters.at("bondlength").get<Real>();
			if (bond_length < 1e-11 || bond_length > 1e-8) {std::cout <<" bondlength out of range 1e-11..1e-8 " << std::endl; success=false;}
		}
		bond_length/=fjc;

		lattice_type=simple_cubic;
		const std::string value = parameters.value("lattice_type", std::string{});
		if (!value.empty()) {
			if (value == "simple_cubic") {lattice_type=simple_cubic; lambda=1.0/6.0; Z=6;}
			else if (value == "hexagonal") {lattice_type=hexagonal; lambda=1.0/4.0; Z=4;}
			else {std::cout << "Input for 'lattice_type' not recognized. 'simple_cubic' or 'hexagonal'." << std::endl; success = false;}
		} else {
			success=false; std::cout <<"Namics can not run without input for 'lattice_type'" << std::endl;
		}

		offset_first_layer =0;
		gradients=parameters.value("gradients",1);
		if (gradients<0||gradients>3) {std::cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << std::endl; gradients=1;}
		success = CheckLatticeInput(parameters) && success;

		if ((fjc>1) && (lattice_type != hexagonal)) {success = false; std::cout << "For FJC-choices >3, we need lattice_type = 'hexagonal'." << std::endl; }
		if (gradients ==2 && fjc>3) {success = false; std::cout <<" When gradients is 2, FJC-choices are limited to 7 " << std::endl; }
		if (gradients ==3 && fjc>2) {success = false; std::cout <<" When gradients is 3, FJC-choices are limited to 5 " << std::endl; }

		if (parameters.contains("ignore_site_fraction")) {
			ignore_sites=parameters.at("ignore_site_fraction").get<bool>();
			if (!ignore_sites) std::cout <<"ignore_site_fraction is set to false. Full site fractions computed. " << std::endl;
		}

		if (parameters.contains("fcc_site_fraction")) {
			fcc_sites=parameters.at("fcc_site_fraction").get<bool>();
			if (!fcc_sites) std::cout <<"fcc_site_fraction is set to false. Full site fractions computed. " << std::endl;
		}

		if (fcc_sites&&ignore_sites) {
			std::cout <<"can't combine 'fcc_site_fraction' with 'ignore_site_fraction'" <<std::endl; success=false;
		}
		stencil_full=true;
		if (parameters.contains("stencil_full")) {
			stencil_full=parameters.at("stencil_full").get<bool>();
			if (gradients<3 && stencil_full) std::cout << "untested territory for 'stencil_full' " << std::endl;
		}
		PutM();
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in lat '" << name << "': " << error.what() << std::endl;
		success = false;
	}

	return success;
}

void Lattice::PrepareForCalculations(void) {
NAMICS_DBG("PrepareForCalculations in lattice" << std::endl);
}

std::span<Real> Lattice::GetPointer(int profile) {
NAMICS_DBG("GetPointer for lattice " + name << std::endl);
	if (profile == 0) return L;
	return {};
}

ParameterStore Lattice::FormatProfile(std::span<const Real> profile, bool write_bounds) {
	const int a = write_bounds ? 0 : fjc;
	ParameterStore out = ParameterStore::array();
	switch (gradients) {
		case 1:
			for (int x = a; x < MX + 2 * fjc - a; ++x) out.push_back(profile[x]);
			break;
		case 2:
			for (int x = a; x < MX + 2 * fjc - a; ++x) {
				ParameterStore row = ParameterStore::array();
				for (int y = a; y < MY + 2 * fjc - a; ++y) row.push_back(profile[P(x, y)]);
				out.push_back(std::move(row));
			}
			break;
		case 3:
			for (int x = a; x < MX + 2 * fjc - a; ++x) {
				ParameterStore plane = ParameterStore::array();
				for (int y = a; y < MY + 2 * fjc - a; ++y) {
					ParameterStore row = ParameterStore::array();
					for (int z = a; z < MZ + 2 * fjc - a; ++z) row.push_back(profile[P(x, y, z)]);
					plane.push_back(std::move(row));
				}
				out.push_back(std::move(plane));
			}
			break;
		default:
			break;
	}
	return out;
}

void Lattice::PushOutput() {
NAMICS_DBG("PushOutput in lat " << std::endl);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["geometry"] = geometry;
	OUTPUT["gradients"] = gradients;
	if (offset_first_layer > 0) OUTPUT["offset_first_layer"] = offset_first_layer;
	OUTPUT["volume"] = volume;
	OUTPUT["accessible volume"] = Accesible_volume;
	OUTPUT["lattice_type"] = lattice_type == simple_cubic ? "simple_cubic" : "hexagonal";
	OUTPUT["bond_length"] = bond_length;
	OUTPUT["FJC_choices"] = FJC;
	OUTPUT["L"] = {{"profile", 0}};

	switch (gradients) {
		case 3:
			OUTPUT["n_layers_z"] = MZ / fjc;

			// Fall through
		case 2:
			OUTPUT["n_layers_y"] = MY / fjc;
			// Fall through
		case 1:
			OUTPUT["n_layers"] = MX / fjc;
			break;
		default:
			break;
	}
}

