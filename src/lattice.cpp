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

std::string GetParameter(const ParameterStore& parameters, const std::string& key) {
	const auto it = parameters.find(key);
	return it == parameters.end() ? std::string() : it->second;
}

bool ParseSelection(const Input& in, const std::string& name, int start, LatticeSelection& selection) {
	ParameterStore parameters;
	if (!in.CheckParameters("lat", name, start, LatticeKeys(), parameters)) return false;

	selection.gradients = ParseInt(GetParameter(parameters, "gradients"), 1);
	if (selection.gradients < 1 || selection.gradients > 3) {
		std::cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << std::endl;
		selection.gradients = 1;
	}

	std::vector<std::string> options = {"spherical", "cylindrical", "flat", "planar"};
	const std::string geometry = GetParameter(parameters, "geometry");
	if (!geometry.empty()) {
		if (!ParseString(geometry, selection.geometry, options, "In lattice input for 'geometry' not recognized."))
			return false;
	} else {
		selection.geometry = "planar";
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
	BC(6), // boundary condition slots: lower/upper for x, y, z
	KEYS(LatticeKeys())
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

	success = In->CheckParameters("lat",name,start, KEYS, PARAMETERS);
	if (!success) return success;
		std::vector<std::string> options;
		FJC=3;	fjc=1;
		if (success && GetValue("FJC_choices").length()>0) {
			if (!ParseInt(GetValue("FJC_choices"),FJC,"FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ..."))
				success=false;
			else if ((FJC-3) %2 != 0) {
				std::cout << "FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ...." <<std::endl;
				success=false;
			}

			fjc=(FJC-1)/2;
		}

		if (success && GetValue("b/l").length()>0) {
			int fjc_new=1;
			if (!ParseInt(GetValue("b/l"),fjc_new,"b/l can adopt only few integer values: 1, 2, 3, ..."))
				success=false;
			else {
				if (fjc_new <1 ) {
					std::cout << "b/l should be a positive integer: 1, 2, 3, ...." <<std::endl;
					success=false;
				}
				if (GetValue("FJC_choices").length()>0 && fjc_new !=fjc) {
					std::cout <<"You have set both 'FJC_choices' and 'b/l', but their values are not consistent with each other."<<std::endl;
					if (fjc_new<fjc && fjc_new >0) {
							std::cout <<"The value of 'b/l' is used, and that of FJC_choices is rejected." << std::endl;
					} else {
							if (fjc_new<1) success=true;
							std::cout <<"The value of 'FJC_choices' is used, and that of b/l is rejected." << std::endl;
							fjc_new=fjc;
					}
				}
			}
			fjc=fjc_new;
			FJC=2*fjc+1; //3+(fjc-1)*2;
		}

		bond_length=0;
		if (GetValue("bondlength").size()>0) {
			bond_length =  ParseReal(GetValue("bondlength"),5e-10);
			if (bond_length < 1e-11 || bond_length > 1e-8) {std::cout <<" bondlength out of range 1e-11..1e-8 " << std::endl; success=false;}
		}
		bond_length/=fjc;

		std::string lat_type;
		lattice_type=simple_cubic;
		options.push_back("simple_cubic"); options.push_back("hexagonal");
		Value=GetValue("lattice_type");
		if (Value.length()>0) {
			if (!ParseString(Value,lat_type,options,"Input for 'lattice_type' not recognized. 'simple_cubic' or 'hexagonal'.")) success = false; else {
				if (lat_type == "simple_cubic") {lattice_type=simple_cubic; lambda=1.0/6.0; Z=6;}
				if (lat_type == "hexagonal") {lattice_type=hexagonal; lambda=1.0/4.0; Z=4;}
			}
		} else {
			success=false; std::cout <<"Namics can not run without input for 'lattice_type'" << std::endl;
		}

		offset_first_layer =0;
		gradients=1;
		gradients=ParseInt(GetValue("gradients"),1);
		if (gradients<0||gradients>3) {std::cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << std::endl; gradients=1;}
		switch(gradients) {
			case 1:
				MX = ParseInt(GetValue("n_layers"),-123);

				if (MX==-123) {success=false; std::cout <<"In 'lat' the parameter 'n_layers' is required. Problem terminated" << std::endl;}
				else {
					if (MX<0 || MX >1e6) {
						success = false;
						std::cout <<"n_layers out of bounds, currently: 0..1e6; Problem terminated" << std::endl;
					}
				}
				MX=fjc*(MX);
				options.clear();
				options.push_back("spherical");
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");

				if (GetValue("geometry").size()>0) {
					if (!ParseString(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized."))
						success=false;
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar";
				if (geometry!="planar") {
					offset_first_layer=ParseReal(GetValue("offset_first_layer"),0);
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
				if (GetValue("lowerbound_x").size()>0) {success=false; std::cout << "lowerbound_x is not allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("lowerbound_y").size()>0) {success=false; std::cout << "lowerbound_y is not allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("lowerbound_z").size()>0) {success=false; std::cout << "lowerbound_z is not allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("upperbound_x").size()>0) {success=false; std::cout << "upperbound_x is not allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("upperbound_y").size()>0) {success=false; std::cout << "upperbound_y is not allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("upperbound_z").size()>0) {success=false; std::cout << "upperbound_z is not allowed in 1-gradient calculations" << std::endl;}


				if (GetValue("lowerbound").size()==0) BC[0]="mirror";
				else if (!ParseString(GetValue("lowerbound"),BC[0],options,"For 'lowerbound' boundary condition not recognized. ")) success=false;

				if (GetValue("upperbound").size()==0) BC[3]="mirror";
				else if (!ParseString(GetValue("upperbound"),BC[3],options,"For 'upperbound' boundary condition not recognized."))
					success = false;

				break;
			case 2:
				if (GetValue("upperbound").size()>0) {success=false; std::cout << "upperbound is only allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("lowerbound").size()>0) {success=false; std::cout << "lowerbound is only allowed in 1-gradient calculations" << std::endl;}

				MX = ParseInt(GetValue("n_layers_x"),-123);
				if (MX==-123) {
					success=false;
					std::cout <<"In 'lat' the parameter 'n_layers_x' is required. Problem terminated" << std::endl;
				}
				else {
					if (MX<0 || MX >1e6) {
						success = false;
						std::cout <<"n_layers_x out of bounds, currently: 0.. 1e6; Problem terminated" << std::endl;
					}
				}
				MX=fjc*(MX);
				MY = ParseInt(GetValue("n_layers_y"),-123);
				if (MY==-123) {
					success=false;
					std::cout <<"In 'lat' the parameter 'n_layers_y' is required. Problem terminated" << std::endl;
				}
				else {
					if (MY<0 || MY >1e6) {
						success = false;
						std::cout <<"n_layers_y out of bounds, currently: 0.. 1e6; Problem terminated" << std::endl;
					}
				}
				MY=fjc*(MY);
				options.clear();
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");
				if (GetValue("geometry").size()>0) {
					if (!ParseString(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized."))
						success=false;
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar";
				if (geometry=="planar") {volume = MX*MY;}

				if (geometry!="planar") {
					offset_first_layer=ParseReal(GetValue("offset_first_layer"),0);
					if (offset_first_layer<0) {
						std::cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << std::endl;
						offset_first_layer=0;
					}
				}

				options.clear();
				options.push_back("mirror");
				options.push_back("surface"); //turned on...let's hope it works
				if (geometry=="planar")
					options.push_back("periodic");
				if (GetValue("lowerbound_z").size()>0) {
					std::cout << "lowerbound_z is not allowed in 2-gradient calculations" << std::endl;
					success=false;
				}
				if (GetValue("upperbound_z").size()>0) {
					std::cout << "upperbound_z is not allowed in 2-gradient calculations" << std::endl;
					success=false;
				}

				if (GetValue("lowerbound_x").size()==0) BC[0]="mirror";
				else if (!ParseString(GetValue("lowerbound_x"),BC[0],options,"for 'lowerbound_x' boundary condition not recognized.  ")) success=false;

				if (GetValue("upperbound_x").size()==0) BC[3]="mirror";
				else if (!ParseString(GetValue("upperbound_x"),BC[3],options,"for 'upperbound_x' boundary condition not recognized. ")) success = false;

				if (GetValue("lowerbound_y").size()==0) BC[1]="mirror";
				else if (!ParseString(GetValue("lowerbound_y"),BC[1],options,"for 'lowerbound_y' boundary condition not recognized. ")) success=false;

				if (GetValue("upperbound_y").size()==0) BC[4]="mirror";
				else if (!ParseString(GetValue("upperbound_y"),BC[4],options,"for 'upperbound_Y' boundary condition not recognized. ")) success = false;


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
			case 3:
				if (GetValue("upperbound").size()>0) {success=false; std::cout << "upperbound is only allowed in 1-gradient calculations" << std::endl;}
				if (GetValue("lowerbound").size()>0) {success=false; std::cout << "lowerbound is only allowed in 1-gradient calculations" << std::endl;}

				if (!ParseInt(GetValue("n_layers_x"),MX,1,1e6,"In 'lat' the parameter 'n_layers_x' is required"))
					success=false;
				if (!ParseInt(GetValue("n_layers_y"),MY,1,1e6,"In 'lat' the parameter 'n_layers_y' is required"))
					success=false;
				if (!ParseInt(GetValue("n_layers_z"),MZ,1,1e6,"In 'lat' the parameter 'n_layers_z' is required"))
					success=false;
				MX=fjc*(MX);
				MY=fjc*(MY);
				MZ=fjc*(MZ);
				options.clear();
				options.push_back("mirror"); //options.push_back("mirror_2");
				options.push_back("periodic");
				//options.push_back("surface");
				//options.push_back("shifted_mirror");


				if (GetValue("lowerbound_x").size()==0) BC[0]="mirror";
				else if (!ParseString(GetValue("lowerbound_x"),BC[0],options,"for 'lowerbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success=false;

				if (GetValue("upperbound_x").size()==0) BC[3]="mirror";
				else if (!ParseString(GetValue("upperbound_x"),BC[3],options,"for 'upperbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success = false;

				if (GetValue("lowerbound_y").size()==0) BC[1]="mirror";
				else if (!ParseString(GetValue("lowerbound_y"),BC[1],options,"for 'lowerbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success=false;

				if (GetValue("upperbound_y").size()==0) BC[4]="mirror";
				else if (!ParseString(GetValue("upperbound_y"),BC[4],options,"for 'upperbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success = false;

				if (GetValue("lowerbound_z").size()==0) BC[2]="mirror";
				else if (!ParseString(GetValue("lowerbound_z"),BC[2],options,"for 'lowerbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success=false;

				if (GetValue("upperbound_z").size()==0) BC[5]="mirror";
				else if (!ParseString(GetValue("upperbound_z"),BC[5],options,"for 'upperbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success = false;

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
			default:
				std::cout << "gradients out of bounds " << std::endl;
				break;
		}


		if ((fjc>1) && (lattice_type != hexagonal)) {success = false; std::cout << "For FJC-choices >3, we need lattice_type = 'hexagonal'." << std::endl; }
		if (gradients ==2 && fjc>3) {success = false; std::cout <<" When gradients is 2, FJC-choices are limited to 7 " << std::endl; }
		if (gradients ==3 && fjc>2) {success = false; std::cout <<" When gradients is 3, FJC-choices are limited to 5 " << std::endl; }

		if (GetValue("ignore_site_fraction").length()>0) {
			ignore_sites=ParseBool(GetValue("ignore_site_fraction"),false);
			if (!ignore_sites) std::cout <<"ignore_site_fraction is set to false. Full site fractions computed. " << std::endl;
		}

		if (GetValue("fcc_site_fraction").length()>0) {
			fcc_sites=ParseBool(GetValue("fcc_site_fraction"),false);
			if (!fcc_sites) std::cout <<"fcc_site_fraction is set to false. Full site fractions computed. " << std::endl;
		}

		if (fcc_sites&&ignore_sites) {
			std::cout <<"can't combine 'fcc_site_fraction' with 'ignore_site_fraction'" <<std::endl; success=false;
		}
		stencil_full=true;
		if (GetValue("stencil_full").length()>0) {
			stencil_full=ParseBool(GetValue("stencil_full"),true);
			if (gradients<3 && stencil_full) std::cout << "untested territory for 'stencil_full' " << std::endl;
		}
		// Initialize system size and indexing.
		PutM();
		Markov=1;

	return success;
}

void Lattice::PutParameter(std::string new_param) {
NAMICS_DBG("PutParameters in lattice " << std::endl); KEYS.push_back(new_param);
}

std::string Lattice::GetValue(std::string parameter){
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

Real Lattice::GetValue(std::span<const Real> X,std::string s){
NAMICS_DBG("GetValue in lattice " << std::endl);if (X.empty()) std::cout << "profile std::span is empty" << std::endl;
	int x=0,y=0,z=0;
	std::vector<std::string> sub;
	In->split(s,',',sub);
	switch(gradients) {
		case 1:
			if (sub.size()==1) {
				x=ParseInt(sub[0],x);
				if (x==-1) x=MX; //trick to get the value of lastlayer; currently only in 1gradient case....
					if (x<0||x>MX+1) {
						std::cout <<"Requested output position is out of bounds." << std::endl;
						return 0;
					} else return X[x];
				} else std::cout <<"Request for profile output does not contain the expected coordinate." << std::endl;
				break;
		case 2:
			if (sub.size()==2) {
				x=ParseInt(sub[0],x);
				y=ParseInt(sub[1],y);
					if (x<0||x>MX+1||y<0||y>MY+1) {
						std::cout <<"Requested output position is out of bounds." << std::endl;
						return 0;
					} else return X[JX*x+y];
				} else std::cout <<"Request for profile output does not contain the expected coordinate." << std::endl;
				break;
		case 3:
			if (sub.size()>2) {
				x=ParseInt(sub[0],x);
				y=ParseInt(sub[1],y);
				z=ParseInt(sub[2],y);
					if (x<0||x>MX+1||y<0||y>MY+1||z<0||z>MZ+1) {
						std::cout <<"Requested output position is out of bounds." << std::endl;
						return 0;
					} else return X[JX*x+JY*y+z];
				} else  std::cout <<"Request for profile output does not contain the expected coordinate." << std::endl;
			return 0;
			break;
		default:
			break;
	}
	return 0;
}



bool Lattice::PrepareForCalculations(void) {
NAMICS_DBG("PrepareForCalculations in lattice" << std::endl);	bool success=true;
	return success;
}

void Lattice::push(std::string s, Real X) {
NAMICS_DBG("push (Real) in lattice " << std::endl); Reals.push_back(s);
	Reals_value.push_back(X);
}
void Lattice::push(std::string s, int X) {
NAMICS_DBG("push (int) in lattice " << std::endl); ints.push_back(s);
	ints_value.push_back(X);
}
void Lattice::push(std::string s, bool X) {
NAMICS_DBG("push (bool) in lattice " << std::endl); bools.push_back(s);
	bools_value.push_back(X);
}
void Lattice::push(std::string s, std::string X) {
NAMICS_DBG("push (std::string) in lattice " << std::endl); strings.push_back(s);
	strings_value.push_back(X);
}

std::span<Real> Lattice::GetPointer(std::string s) {
NAMICS_DBG("GetPointer for lattice " + name << std::endl);	std::vector<std::string> sub;
	In->split(s,';',sub);
	if (sub[0]=="profile" && sub[1]=="0") return L;
	if (sub[0]=="std::vector") {}
	return {};
}

std::span<int> Lattice::GetPointerInt(std::string s) {
NAMICS_DBG("GetPointerInt for lattice " + name << std::endl);	std::vector<std::string> sub;
	In->split(s,';',sub);
	if (sub[0]=="std::array"){//get with sub[1] the number and put the pointer to integer std::array in return.
	}

	return {};
}

void Lattice::PushOutput() {
NAMICS_DBG("PushOutput in lat " << std::endl); strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	std::string mirror="mirror";
	std::string periodic="periodic";
	push("geometry",geometry);
	push("gradients",gradients);
	if (offset_first_layer>0) push("offset_first_layer",offset_first_layer);
	push("volume",volume);
	push("accessible volume",Accesible_volume);
	std::string LatticeType;
	if (lattice_type == simple_cubic) LatticeType="simple_cubic"; else LatticeType="hexagonal";
	push("lattice_type",LatticeType);
	push("bond_length",bond_length);
	push("FJC_choices",FJC);
	std::string s="profile;0"; push("L",s);

	switch (gradients) {
		case 3:
			push("n_layers_z",MZ/fjc);

			// Fall through
		case 2:
			push("n_layers_y",MY/fjc);
			// Fall through
		case 1:
			push("n_layers",MX/fjc);
			break;
		default:
			break;
	}
}


int Lattice::GetValue(std::string prop,int &int_result,Real &Real_result,std::string &string_result){
NAMICS_DBG("GetValue (long)  in lattice " << std::endl);
	for ( size_t i = 0 ; i<ints.size() ; ++i)
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}

	for ( size_t i = 0 ; i<Reals.size() ; ++i)
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}

	for ( size_t i = 0 ; i<bools.size() ; ++i)
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true";
			else string_result="false";
			return 3;
		}

	for ( size_t i = 0 ; i<strings.size() ; ++i)
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}

	return 0;
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
