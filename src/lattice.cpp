#include "lattice.h"
#include "LG1Planar.h"
#include "LG2Planar.h"
#include "LGrad1.h"
#include "LGrad2.h"
#include "LGrad3.h"

namespace {

const std::vector<std::string>& LatticeKeys() {
	static const std::vector<std::string> keys = {
		"gradients", "n_layers", "offset_first_layer", "geometry",
		"n_layers_x", "n_layers_y", "n_layers_z",
		"lowerbound", "upperbound",
		"lowerbound_x", "upperbound_x",
		"lowerbound_y", "upperbound_y",
		"lowerbound_z", "upperbound_z",
		"bondlength", "ignore_site_fraction", "fcc_site_fraction",
		"lattice_type", "stencil_full", "FJC_choices", "b/l"
	};
	return keys;
}

bool ParseSelection(const Input& in, const std::string& name, int start, LatticeSelection& selection) {
	const auto& parameters = in.Parameters("lat", name, start);
	bool success = true;
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(LatticeKeys(), it.key())) continue;
		success = false;
		std::cout << "lat property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : LatticeKeys()) std::cout << item << std::endl;
	}
	if (!success) return false;

	try {
		selection.gradients = parameters.value("gradients", 1);
		if (selection.gradients < 1 || selection.gradients > 3) {
			std::cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << std::endl;
			selection.gradients = 1;
		}

		selection.geometry = parameters.value("geometry", std::string{"planar"});
		if (selection.geometry != "spherical" && selection.geometry != "cylindrical" && selection.geometry != "flat" && selection.geometry != "planar") {
			std::cout << "In lattice input for 'geometry' not recognized." << std::endl;
			return false;
		}
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in lat '" << name << "': " << error.what() << std::endl;
		return false;
	}
	if (selection.geometry == "flat") selection.geometry = "planar";
	return true;
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
	Markov=1;
	subl=0;
}

namespace lattice_factory {

std::unique_ptr<Lattice> CreateChecked(const Input& in, const std::string& name, int start) {
	LatticeSelection selection;
	if (!ParseSelection(in, name, start, selection)) return nullptr;

	if (auto lattice = TryCreate<LG1Planar>(in, name, start, selection)) return lattice;
	if (auto lattice = TryCreate<LGrad1>(in, name, start, selection)) return lattice;
	if (auto lattice = TryCreate<LG2Planar>(in, name, start, selection)) return lattice;
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
	l1.clear();
	l11.clear();
	l_1.clear();
	l_11.clear();
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
	LABDA.clear();
	LABDA_1.clear();
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
			if (BC[2]=="mirror") {
				if (fjc==1) BZ1=1; else {
					for (int k=0; k<fjc; k++) B_Z1[k]=2*fjc-1-k;
				}
			}
			if (BC[2]=="periodic") {
				if (fjc==1) BZ1=MZ; else {
					for (int k=0; k<fjc; k++) B_Z1[k]=MZ+k;
				}
			}
			if (BC[2]=="surface") {
				if (fjc==1) BZ1=0; else {
					for (int k=0; k<fjc; k++) B_Z1[k]=k;
				}
			}

			if (BC[5]=="mirror") {
				if (fjc==1) BZM=MZ; else {
					for (int k=0; k<fjc; k++) B_ZM[k]=MZ+fjc-k-1;
				}
			}
			if (BC[5]=="periodic") {
				if (fjc==1) BZM=1; else {
					for (int k=0; k<fjc; k++) B_ZM[k]=fjc+k;
				}
			}
			if (BC[5]=="surface") {
				if (fjc==1) BZM=MZ+1; else {
					for (int k=0; k<fjc; k++) B_ZM[k]=MZ+fjc+k;
				}
			}

			//Fall through
		case 2:
			if (BC[1]=="mirror") {
				if (fjc==1) BY1=1; else {
					for (int k=0; k<fjc; k++) B_Y1[k]=2*fjc-1-k;
				}
			}
			if (BC[1]=="periodic") {
				if (fjc==1) BY1=MY; else {
					for (int k=0; k<fjc; k++) B_Y1[k]=MY+k;
				}
			}
			if (BC[1]=="surface") {
				if (fjc==1) BY1=0; else {
					for (int k=0; k<fjc; k++) B_Y1[k]=k;
				}
			}

			if (BC[4]=="mirror") {
				if (fjc==1) BYM=MY; else {
					for (int k=0; k<fjc; k++) B_YM[k]=MY+fjc-k-1;
				}
			}
			if (BC[4]=="periodic") {
				if (fjc==1) BYM=1; else {
					for (int k=0; k<fjc; k++) B_YM[k]=fjc+k;
				}
			}
			if (BC[4]=="surface") {
				if (fjc==1) BYM=MY+1; else { //std::cout <<"surface ub" << std::endl;
					for (int k=0; k<fjc; k++) B_YM[k]=MY+fjc+k;
				}
			}

			//Fall through
		case 1:
			if (BC[0]=="mirror") {
				if (fjc==1) BX1=1; else {
					for (int k=0; k<fjc; k++) B_X1[k]=2*fjc-1-k;
				}
			}
			if (BC[0]=="periodic") {
				if (fjc==1) BX1=MX;else {
					for (int k=0; k<fjc; k++) B_X1[k]=MX+k;
				}
			}
			if (BC[0]=="surface") {
				if (fjc==1) BX1=0; else {
					for (int k=0; k<fjc; k++) B_X1[k]=k;
				}
			}

			if (BC[3]=="mirror") {
				if (fjc==1) BXM=MX;else {
					for (int k=0; k<fjc; k++) B_XM[k]=MX+fjc-k-1;
				}
			}
			if (BC[3]=="periodic") {
				if (fjc==1) BXM=1; else {
					for (int k=0; k<fjc; k++) B_XM[k]=fjc+k;
				}
			}
			if (BC[3]=="surface") {
				if (fjc==1) BXM=MX+1; else {
					for (int k=0; k<fjc; k++) B_XM[k]=MX+fjc+k;
				}
			}

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
	if (Markov==2) {
		if (fjc==1) {
			l1.assign(M, 0);
			l_1.assign(M, 0);
			l11.assign(M, 0);
			l_11.assign(M, 0);
		} else {
			LABDA.assign(FJC * M, 0);
			LABDA_1.assign(FJC * M, 0);
		}
		H.assign(M, 0);
	}


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


bool Lattice::PutSub_box(int mx_, int my_, int mz_,int n_box_) {
	bool success = true;
	if (mx_<1 || my_<1 || mz_<1 || mx_>MX || my_>MY || mz_>MZ) {std::cout <<"subbox size out of bound: mx= " << mx_ << " my = " << my_ << " mz = " << mz_ << ", while MX = " << MX << " MY = " << MY << " MZ = " << MZ  << std::endl; success=false; }
	mx.push_back(mx_); my.push_back(my_); mz.push_back(mz_);
	m.push_back((mx_+2)*(my_+2)*(mz_+2));
	jx.push_back((mx_+2)*(my_+2)); jy.push_back(my_+2);
	n_box.push_back(n_box_);
	return success;
}

bool Lattice::CheckInput(int start) {
NAMICS_DBG("CheckInput in lattice " << std::endl);	bool success=true;
	mx.push_back(0); my.push_back(0); mz.push_back(0); jx.push_back(0); jy.push_back(0); m.push_back(0); n_box.push_back(0);
	std::string Value;
	const auto& parameters = In->Parameters("lat", name, start);

	try {
		auto assign_choice = [&](const std::string& value, std::string& target, const std::vector<std::string>& allowed, const std::string& error) {
			if (!ContainsValue(allowed, value)) {
				std::cout << error << std::endl;
				success = false;
				return;
			}
			target = value;
		};
		std::vector<std::string> options;
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

		std::string lat_type;
		lattice_type=simple_cubic;
		Value=parameters.value("lattice_type", std::string{});
		if (Value.length()>0) {
			if (Value == "simple_cubic") {lat_type = Value; lattice_type=simple_cubic; lambda=1.0/6.0; Z=6;}
			else if (Value == "hexagonal") {lat_type = Value; lattice_type=hexagonal; lambda=1.0/4.0; Z=4;}
			else {std::cout << "Input for 'lattice_type' not recognized. 'simple_cubic' or 'hexagonal'." << std::endl; success = false;}
		} else {
			success=false; std::cout <<"Namics can not run without input for 'lattice_type'" << std::endl;
		}

		offset_first_layer =0;
		gradients=parameters.value("gradients",1);
		if (gradients<0||gradients>3) {std::cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << std::endl; gradients=1;}
		switch(gradients) {
			case 1: {
				MX = parameters.contains("n_layers") ? parameters.at("n_layers").get<int>() : -123;
				if (MX==-123) {success=false; std::cout <<"In 'lat' the parameter 'n_layers' is required. Problem terminated" << std::endl;}
				else if (MX<0 || MX >1e6) {
					success = false;
					std::cout <<"n_layers out of bounds, currently: 0..1e6; Problem terminated" << std::endl;
				}
				MX=fjc*(MX);
				options.clear();
				options.push_back("spherical");
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");

				const std::string geometry_value = parameters.value("geometry", std::string{});
				if (!geometry_value.empty()) {
					assign_choice(geometry_value, geometry, options, "In lattice input for 'geometry' not recognized.");
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar";
				if (geometry!="planar") {
					offset_first_layer=parameters.value("offset_first_layer",0.0);
					if (offset_first_layer<0) {
						std::cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << std::endl;
						offset_first_layer=0;
					}
					volume = MX/fjc;
				}
				offset_first_layer *=fjc;

				options.clear();
				options.push_back("mirror");
				options.push_back("surface");
				options.push_back("periodic");
				if (parameters.contains("lowerbound_x")) {success=false; std::cout << "lowerbound_x is not allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("lowerbound_y")) {success=false; std::cout << "lowerbound_y is not allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("lowerbound_z")) {success=false; std::cout << "lowerbound_z is not allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("upperbound_x")) {success=false; std::cout << "upperbound_x is not allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("upperbound_y")) {success=false; std::cout << "upperbound_y is not allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("upperbound_z")) {success=false; std::cout << "upperbound_z is not allowed in 1-gradient calculations" << std::endl;}

				const std::string lowerbound = parameters.value("lowerbound", std::string{});
				if (lowerbound.empty()) BC[0]="mirror";
				else assign_choice(lowerbound, BC[0], options, "For 'lowerbound' boundary condition not recognized. ");

				const std::string upperbound = parameters.value("upperbound", std::string{});
				if (upperbound.empty()) BC[3]="mirror";
				else assign_choice(upperbound, BC[3], options, "For 'upperbound' boundary condition not recognized.");
				break;
			}
			case 2: {
				if (parameters.contains("upperbound")) {success=false; std::cout << "upperbound is only allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("lowerbound")) {success=false; std::cout << "lowerbound is only allowed in 1-gradient calculations" << std::endl;}

				MX = parameters.contains("n_layers_x") ? parameters.at("n_layers_x").get<int>() : -123;
				if (MX==-123) {
					success=false;
					std::cout <<"In 'lat' the parameter 'n_layers_x' is required. Problem terminated" << std::endl;
				} else if (MX<0 || MX >1e6) {
					success = false;
					std::cout <<"n_layers_x out of bounds, currently: 0.. 1e6; Problem terminated" << std::endl;
				}
				MX=fjc*(MX);
				MY = parameters.contains("n_layers_y") ? parameters.at("n_layers_y").get<int>() : -123;
				if (MY==-123) {
					success=false;
					std::cout <<"In 'lat' the parameter 'n_layers_y' is required. Problem terminated" << std::endl;
				} else if (MY<0 || MY >1e6) {
					success = false;
					std::cout <<"n_layers_y out of bounds, currently: 0.. 1e6; Problem terminated" << std::endl;
				}
				MY=fjc*(MY);
				options.clear();
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");
				const std::string geometry_value = parameters.value("geometry", std::string{});
				if (!geometry_value.empty()) {
					assign_choice(geometry_value, geometry, options, "In lattice input for 'geometry' not recognized.");
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar";
				if (geometry=="planar") {volume = MX*MY;}

				if (geometry!="planar") {
					offset_first_layer=parameters.value("offset_first_layer",0.0);
					if (offset_first_layer<0) {
						std::cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << std::endl;
						offset_first_layer=0;
					}
				}

				options.clear();
				options.push_back("mirror");
				options.push_back("surface");
				if (geometry=="planar") options.push_back("periodic");
				if (parameters.contains("lowerbound_z")) {
					std::cout << "lowerbound_z is not allowed in 2-gradient calculations" << std::endl;
					success=false;
				}
				if (parameters.contains("upperbound_z")) {
					std::cout << "upperbound_z is not allowed in 2-gradient calculations" << std::endl;
					success=false;
				}

				const std::string lowerbound_x = parameters.value("lowerbound_x", std::string{});
				if (lowerbound_x.empty()) BC[0]="mirror";
				else assign_choice(lowerbound_x, BC[0], options, "for 'lowerbound_x' boundary condition not recognized.  ");

				const std::string upperbound_x = parameters.value("upperbound_x", std::string{});
				if (upperbound_x.empty()) BC[3]="mirror";
				else assign_choice(upperbound_x, BC[3], options, "for 'upperbound_x' boundary condition not recognized. ");

				const std::string lowerbound_y = parameters.value("lowerbound_y", std::string{});
				if (lowerbound_y.empty()) BC[1]="mirror";
				else assign_choice(lowerbound_y, BC[1], options, "for 'lowerbound_y' boundary condition not recognized. ");

				const std::string upperbound_y = parameters.value("upperbound_y", std::string{});
				if (upperbound_y.empty()) BC[4]="mirror";
				else assign_choice(upperbound_y, BC[4], options, "for 'upperbound_Y' boundary condition not recognized. ");

				if (BC[0]=="periodic" || BC[3]=="periodic") {
					if (BC[0]!=BC[3]) {
						success=false;
						std::cout <<"For boundaries in x-direction: 'periodic' BC  should be set to upper and lower bounds " << std::endl;
					}
				}
				if (BC[1]=="periodic" || BC[4]=="periodic") {
					if (BC[1]!=BC[4]) {success=false;  std::cout <<"For boundaries in y-direction: 'periodic' BC should be set to upper and lower bounds " << std::endl;}
				}
				break;
			}
			case 3: {
				if (parameters.contains("upperbound")) {success=false; std::cout << "upperbound is only allowed in 1-gradient calculations" << std::endl;}
				if (parameters.contains("lowerbound")) {success=false; std::cout << "lowerbound is only allowed in 1-gradient calculations" << std::endl;}

				if (!parameters.contains("n_layers_x")) {std::cout <<"In 'lat' the parameter 'n_layers_x' is required" << std::endl; success=false;}
				else MX = parameters.at("n_layers_x").get<int>();
				if (!parameters.contains("n_layers_y")) {std::cout <<"In 'lat' the parameter 'n_layers_y' is required" << std::endl; success=false;}
				else MY = parameters.at("n_layers_y").get<int>();
				if (!parameters.contains("n_layers_z")) {std::cout <<"In 'lat' the parameter 'n_layers_z' is required" << std::endl; success=false;}
				else MZ = parameters.at("n_layers_z").get<int>();
				if (MX<1 || MX>1e6) success=false;
				if (MY<1 || MY>1e6) success=false;
				if (MZ<1 || MZ>1e6) success=false;
				MX=fjc*(MX);
				MY=fjc*(MY);
				MZ=fjc*(MZ);
				options.clear();
				options.push_back("mirror");
				options.push_back("periodic");

				const std::string lowerbound_x = parameters.value("lowerbound_x", std::string{});
				if (lowerbound_x.empty()) BC[0]="mirror";
				else assign_choice(lowerbound_x, BC[0], options, "for 'lowerbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ");

				const std::string upperbound_x = parameters.value("upperbound_x", std::string{});
				if (upperbound_x.empty()) BC[3]="mirror";
				else assign_choice(upperbound_x, BC[3], options, "for 'upperbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ");

				const std::string lowerbound_y = parameters.value("lowerbound_y", std::string{});
				if (lowerbound_y.empty()) BC[1]="mirror";
				else assign_choice(lowerbound_y, BC[1], options, "for 'lowerbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ");

				const std::string upperbound_y = parameters.value("upperbound_y", std::string{});
				if (upperbound_y.empty()) BC[4]="mirror";
				else assign_choice(upperbound_y, BC[4], options, "for 'upperbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ");

				const std::string lowerbound_z = parameters.value("lowerbound_z", std::string{});
				if (lowerbound_z.empty()) BC[2]="mirror";
				else assign_choice(lowerbound_z, BC[2], options, "for 'lowerbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ");

				const std::string upperbound_z = parameters.value("upperbound_z", std::string{});
				if (upperbound_z.empty()) BC[5]="mirror";
				else assign_choice(upperbound_z, BC[5], options, "for 'upperbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ");

				if (BC[1]=="periodic" || BC[4]=="periodic") {
					if (BC[1] != BC[4]) {
						std::cout <<"In y-direction the boundary conditions do not match:" + BC[1] << " and " <<  BC[4] << std::endl;
						success=false;
					}
				}
				if (BC[2]=="periodic" || BC[5]=="periodic") {
					if (BC[2] != BC[5]) {
						std::cout <<"In z-direction the boundary conditions do not match:" + BC[2] << " and " <<  BC[5] << std::endl;
						success=false;
					}
				}
				break;
			}
			default:
				std::cout << "gradients out of bounds " << std::endl;
				break;
		}

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
		Markov=1;
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in lat '" << name << "': " << error.what() << std::endl;
		success = false;
	}

	return success;
}

bool Lattice::PrepareForCalculations(void) {
NAMICS_DBG("PrepareForCalculations in lattice" << std::endl);	bool success=true;
	return success;
}

std::span<Real> Lattice::GetPointer(int profile) {
NAMICS_DBG("GetPointer for lattice " + name << std::endl);
	if (profile == 0) return L;
	return {};
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




void Lattice::DistributeG1(std::span<const Real> G1, std::span<Real> g1, std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int n_box) {
	int k=sub_box_on;
	tools::DistributeG1(G1, g1, Bx, By, Bz, M, m[k], n_box, mx[k], my[k], mz[k], MX, MY, MZ, jx[k], jy[k], JX, JY);
}

void Lattice::CollectPhi(std::span<Real> phi, std::span<const Real> GN, std::span<const Real> rho, std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int n_box) {
	int k=sub_box_on;
	tools::CollectPhi(phi, GN, rho, Bx, By, Bz, M, m[k], n_box, mx[k], my[k], mz[k], MX, MY, MZ, jx[k], jy[k], JX, JY);
}

void Lattice::ComputeGN(std::span<Real> GN, std::span<const Real> Gg_f, std::span<const int> H_Bx, std::span<const int> H_By, std::span<const int> H_Bz, std::span<const int> H_Px2, std::span<const int> H_Py2, std::span<const int> H_Pz2, int N, int n_box) {
	int k=sub_box_on;
	for (int p=0; p<n_box; p++) std::copy_n(Gg_f.data()+n_box*m[k]*N +p*m[k]+ jx[k]*(H_Px2[p]-H_Bx[p])+jy[k]*(H_Py2[p]-H_By[p])+(H_Pz2[p]-H_Bz[p]), 1, GN.data()+p);

}
