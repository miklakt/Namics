#include "lattice.h"

Lattice::Lattice(const Input& In_,const string& name_) :
	BC(6) // boundary condition slots: lower/upper for x, y, z
{ //this file contains switch (gradients). In this way we keep all the lattice issues in one file!
NAMICS_DBG("Lattice constructor" << endl);	In=&In_; name=name_;
	KEYS.push_back("gradients"); KEYS.push_back("n_layers"); KEYS.push_back("offset_first_layer");
	KEYS.push_back("geometry");
	KEYS.push_back("n_layers_x");   KEYS.push_back("n_layers_y"); KEYS.push_back("n_layers_z");
	KEYS.push_back("lowerbound"); KEYS.push_back("upperbound");
	KEYS.push_back("lowerbound_x"); KEYS.push_back("upperbound_x");
	KEYS.push_back("lowerbound_y"); KEYS.push_back("upperbound_y");
	KEYS.push_back("lowerbound_z"); KEYS.push_back("upperbound_z");
 	KEYS.push_back("bondlength");
	KEYS.push_back("ignore_site_fraction");
	KEYS.push_back("fcc_site_fraction");
  	KEYS.push_back("lattice_type");
	KEYS.push_back("stencil_full");
	KEYS.push_back("FJC_choices");
	KEYS.push_back("b/l");
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
	l1 = nullptr;
	l11 = nullptr;
	l_1 = nullptr;
	l_11 = nullptr;
	H = nullptr;
	B_X1 = nullptr;
	B_Y1 = nullptr;
	B_Z1 = nullptr;
	B_XM = nullptr;
	B_YM = nullptr;
	B_ZM = nullptr;
	L = nullptr;
	lambda0 = nullptr;
	fcc_lambda0 = nullptr;
	lambda_1 = nullptr;
	fcc_lambda_1 = nullptr;
	lambda1 = nullptr;
	fcc_lambda1 = nullptr;
	LAMBDA = nullptr;
	LABDA = nullptr;
	LABDA_1 = nullptr;
	X = nullptr;
}

void Lattice::DeAllocateMemory(void) {
NAMICS_DBG("DeAllocateMemory in lat " << endl);	if (!all_lattice) return;
	all_lattice=false;
	delete[] l1;
	delete[] l11;
	delete[] l_1;
	delete[] l_11;
	delete[] H;
	delete[] B_X1;
	delete[] B_Y1;
	delete[] B_Z1;
	delete[] B_XM;
	delete[] B_YM;
	delete[] B_ZM;
	delete[] L;
	delete[] lambda0;
	delete[] fcc_lambda0;
	delete[] lambda_1;
	delete[] fcc_lambda_1;
	delete[] lambda1;
	delete[] fcc_lambda1;
	delete[] LAMBDA;
	delete[] LABDA;
	delete[] LABDA_1;
	delete[] X;
	l1 = nullptr;
	l11 = nullptr;
	l_1 = nullptr;
	l_11 = nullptr;
	H = nullptr;
	B_X1 = nullptr;
	B_Y1 = nullptr;
	B_Z1 = nullptr;
	B_XM = nullptr;
	B_YM = nullptr;
	B_ZM = nullptr;
	L = nullptr;
	lambda0 = nullptr;
	fcc_lambda0 = nullptr;
	lambda_1 = nullptr;
	fcc_lambda_1 = nullptr;
	lambda1 = nullptr;
	fcc_lambda1 = nullptr;
	LAMBDA = nullptr;
	LABDA = nullptr;
	LABDA_1 = nullptr;
	X = nullptr;
}


void Lattice::AllocateMemory(void) {
NAMICS_DBG("AllocateMemory in lat " << endl);
	DeAllocateMemory();
	all_lattice=true;
	PutM();
	if (fjc>1) {
		B_X1 = new int[fjc]();
		B_Y1 = new int[fjc]();
		B_Z1 = new int[fjc]();
		B_XM = new int[fjc]();
		B_YM = new int[fjc]();
		B_ZM = new int[fjc]();
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
				if (fjc==1) BYM=MY+1; else { //cout <<"surface ub" << endl;
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
		fcc_lambda_1 = new Real[M]();
		fcc_lambda1 = new Real[M]();
		fcc_lambda0 = new Real[M]();
	}

	if (fjc==1) {
		if (gradients<3) {
		L = new Real[M]();
		lambda_1 = new Real[M]();
		lambda1 = new Real[M]();
		lambda0 = new Real[M]();
		}
	} else {
		L = new Real[M]();
		LAMBDA = new Real[FJC*M]();
	}
	if (Markov==2) {
		if (fjc==1) {
			l1 = new Real[M]();
			l_1 = new Real[M]();
			l11 = new Real[M]();
			l_11 = new Real[M]();
		} else {
			LABDA = new Real[FJC*M]();
			LABDA_1 = new Real[FJC*M]();
		}
		H = new Real[M]();
	}


	X = new Real[M]();
	ComputeLambdas();
}



Lattice::~Lattice() {
NAMICS_DBG("lattice destructor " << endl); DeAllocateMemory();

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
	if (mx_<1 || my_<1 || mz_<1 || mx_>MX || my_>MY || mz_>MZ) {cout <<"subbox size out of bound: mx= " << mx_ << " my = " << my_ << " mz = " << mz_ << ", while MX = " << MX << " MY = " << MY << " MZ = " << MZ  << endl; success=false; }
	mx.push_back(mx_); my.push_back(my_); mz.push_back(mz_);
	m.push_back((mx_+2)*(my_+2)*(mz_+2));
	jx.push_back((mx_+2)*(my_+2)); jy.push_back(my_+2);
	n_box.push_back(n_box_);
	return success;
}

bool Lattice::CheckInput(int start, bool checking) {
NAMICS_DBG("CheckInput in lattice " << endl);	bool success=true;
	mx.push_back(0); my.push_back(0); mz.push_back(0); jx.push_back(0); jy.push_back(0); m.push_back(0); n_box.push_back(0);
	string Value;

	success = In->CheckParameters("lat",name,start, KEYS, PARAMETERS);
	if (!success) return success;
		vector<string> options;
		if (checking) {
			gradients=1;
			gradients=ParseInt(GetValue("gradients"),1);
			if (gradients<0||gradients>3) {cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << endl; gradients=1;}
			options.clear();
			options.push_back("spherical");
			options.push_back("cylindrical");
			options.push_back("flat");options.push_back("planar");

			if (GetValue("geometry").size()>0) {
				if (!ParseString(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized."))
					success=false;
			} else geometry = "planar";
			if (geometry=="flat") geometry="planar";

			return success;
		}



		FJC=3;	fjc=1;
		if (success && GetValue("FJC_choices").length()>0) {
			if (!ParseInt(GetValue("FJC_choices"),FJC,"FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ..."))
				success=false;
			else if ((FJC-3) %2 != 0) {
				cout << "FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ...." <<endl;
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
					cout << "b/l should be a positive integer: 1, 2, 3, ...." <<endl;
					success=false;
				}
				if (GetValue("FJC_choices").length()>0 && fjc_new !=fjc) {
					cout <<"You have set both 'FJC_choices' and 'b/l', but their values are not consistent with each other."<<endl;
					if (fjc_new<fjc && fjc_new >0) {
							cout <<"The value of 'b/l' is used, and that of FJC_choices is rejected." << endl;
					} else {
							if (fjc_new<1) success=true;
							cout <<"The value of 'FJC_choices' is used, and that of b/l is rejected." << endl;
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
			if (bond_length < 1e-11 || bond_length > 1e-8) {cout <<" bondlength out of range 1e-11..1e-8 " << endl; success=false;}
		}
		bond_length/=fjc;

		string lat_type;
		lattice_type=simple_cubic;
		options.push_back("simple_cubic"); options.push_back("hexagonal");
		Value=GetValue("lattice_type");
		if (Value.length()>0) {
			if (!ParseString(Value,lat_type,options,"Input for 'lattice_type' not recognized. 'simple_cubic' or 'hexagonal'.")) success = false; else {
				if (lat_type == "simple_cubic") {lattice_type=simple_cubic; lambda=1.0/6.0; Z=6;}
				if (lat_type == "hexagonal") {lattice_type=hexagonal; lambda=1.0/4.0; Z=4;}
			}
		} else {
			success=false; cout <<"Namics can not run without input for 'lattice_type'" << endl;
		}

		offset_first_layer =0;
		gradients=1;
		gradients=ParseInt(GetValue("gradients"),1);
		if (gradients<0||gradients>3) {cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << endl; gradients=1;}
		switch(gradients) {
			case 1:
				MX = ParseInt(GetValue("n_layers"),-123);

				if (MX==-123) {success=false; cout <<"In 'lat' the parameter 'n_layers' is required. Problem terminated" << endl;}
				else {
					if (MX<0 || MX >1e6) {
						success = false;
						cout <<"n_layers out of bounds, currently: 0..1e6; Problem terminated" << endl;
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
						cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << endl;
						offset_first_layer=0;
					}
					volume = MX/fjc;
				}
				offset_first_layer *=fjc;

				options.clear();
				options.push_back("mirror");
				options.push_back("surface");
				options.push_back("periodic");
				if (GetValue("lowerbound_x").size()>0) {success=false; cout << "lowerbound_x is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound_y").size()>0) {success=false; cout << "lowerbound_y is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound_z").size()>0) {success=false; cout << "lowerbound_z is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_x").size()>0) {success=false; cout << "upperbound_x is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_y").size()>0) {success=false; cout << "upperbound_y is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_z").size()>0) {success=false; cout << "upperbound_z is not allowed in 1-gradient calculations" << endl;}


				if (GetValue("lowerbound").size()==0) BC[0]="mirror";
				else if (!ParseString(GetValue("lowerbound"),BC[0],options,"For 'lowerbound' boundary condition not recognized. ")) success=false;

				if (GetValue("upperbound").size()==0) BC[3]="mirror";
				else if (!ParseString(GetValue("upperbound"),BC[3],options,"For 'upperbound' boundary condition not recognized."))
					success = false;

				break;
			case 2:
				if (GetValue("upperbound").size()>0) {success=false; cout << "upperbound is only allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound").size()>0) {success=false; cout << "lowerbound is only allowed in 1-gradient calculations" << endl;}

				MX = ParseInt(GetValue("n_layers_x"),-123);
				if (MX==-123) {
					success=false;
					cout <<"In 'lat' the parameter 'n_layers_x' is required. Problem terminated" << endl;
				}
				else {
					if (MX<0 || MX >1e6) {
						success = false;
						cout <<"n_layers_x out of bounds, currently: 0.. 1e6; Problem terminated" << endl;
					}
				}
				MX=fjc*(MX);
				MY = ParseInt(GetValue("n_layers_y"),-123);
				if (MY==-123) {
					success=false;
					cout <<"In 'lat' the parameter 'n_layers_y' is required. Problem terminated" << endl;
				}
				else {
					if (MY<0 || MY >1e6) {
						success = false;
						cout <<"n_layers_y out of bounds, currently: 0.. 1e6; Problem terminated" << endl;
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
						cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << endl;
						offset_first_layer=0;
					}
				}

				options.clear();
				options.push_back("mirror");
				options.push_back("surface"); //turned on...let's hope it works
				if (geometry=="planar")
					options.push_back("periodic");
				if (GetValue("lowerbound_z").size()>0) {
					cout << "lowerbound_z is not allowed in 2-gradient calculations" << endl;
					success=false;
				}
				if (GetValue("upperbound_z").size()>0) {
					cout << "upperbound_z is not allowed in 2-gradient calculations" << endl;
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
						cout <<"For boundaries in x-direction: 'periodic' BC  should be set to upper and lower bounds " << endl;
					}
				}
				if (BC[1]=="periodic" || BC[4]=="periodic") {
					if (BC[1]!=BC[4]) {success=false;  cout <<"For boundaries in y-direction: 'periodic' BC should be set to upper and lower bounds " << endl;}
				}
				break;
			case 3:
				if (GetValue("upperbound").size()>0) {success=false; cout << "upperbound is only allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound").size()>0) {success=false; cout << "lowerbound is only allowed in 1-gradient calculations" << endl;}

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
					cout <<"In y-direction the boundary conditions do not match:" + BC[1] << " and " <<  BC[4] << endl;
					success=false;
					}
				}
				if (BC[2]=="periodic" || BC[5]=="periodic") {
					if (BC[2] != BC[5]) {
						cout <<"In z-direction the boundary conditions do not match:" + BC[2] << " and " <<  BC[5] << endl;
						success=false;
					}
				}

				break;
			default:
				cout << "gradients out of bounds " << endl;
				break;
		}


		if ((fjc>1) && (lattice_type != hexagonal)) {success = false; cout << "For FJC-choices >3, we need lattice_type = 'hexagonal'." << endl; }
		if (gradients ==2 && fjc>3) {success = false; cout <<" When gradients is 2, FJC-choices are limited to 7 " << endl; }
		if (gradients ==3 && fjc>2) {success = false; cout <<" When gradients is 3, FJC-choices are limited to 5 " << endl; }

		if (GetValue("ignore_site_fraction").length()>0) {
			ignore_sites=ParseBool(GetValue("ignore_site_fraction"),false);
			if (!ignore_sites) cout <<"ignore_site_fraction is set to false. Full site fractions computed. " << endl;
		}

		if (GetValue("fcc_site_fraction").length()>0) {
			fcc_sites=ParseBool(GetValue("fcc_site_fraction"),false);
			if (!fcc_sites) cout <<"fcc_site_fraction is set to false. Full site fractions computed. " << endl;
		}

		if (fcc_sites&&ignore_sites) {
			cout <<"can't combine 'fcc_site_fraction' with 'ignore_site_fraction'" <<endl; success=false;
		}
		stencil_full=true;
		if (GetValue("stencil_full").length()>0) {
			stencil_full=ParseBool(GetValue("stencil_full"),true);
			if (gradients<3 && stencil_full) cout << "untested territory for 'stencil_full' " << endl;
		}
		// Initialize system size and indexing.
		PutM();
		Markov=1;

	return success;
}

void Lattice::PutParameter(string new_param) {
NAMICS_DBG("PutParameters in lattice " << endl); KEYS.push_back(new_param);
}

string Lattice::GetValue(string parameter){
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

Real Lattice::GetValue(Real* X,string s){
NAMICS_DBG("GetValue in lattice " << endl);if (X==NULL) cout << "pointer X is zero" << endl;
	int x=0,y=0,z=0;
	vector<string> sub;
	In->split(s,',',sub);
	switch(gradients) {
		case 1:
			if (sub.size()==1) {
				x=ParseInt(sub[0],x);
				if (x==-1) x=MX; //trick to get the value of lastlayer; currently only in 1gradient case....
					if (x<0||x>MX+1) {
						cout <<"Requested output position is out of bounds." << endl;
						return 0;
					} else return X[x];
				} else cout <<"Request for profile output does not contain the expected coordinate." << endl;
				break;
		case 2:
			if (sub.size()==2) {
				x=ParseInt(sub[0],x);
				y=ParseInt(sub[1],y);
					if (x<0||x>MX+1||y<0||y>MY+1) {
						cout <<"Requested output position is out of bounds." << endl;
						return 0;
					} else return X[JX*x+y];
				} else cout <<"Request for profile output does not contain the expected coordinate." << endl;
				break;
		case 3:
			if (sub.size()>2) {
				x=ParseInt(sub[0],x);
				y=ParseInt(sub[1],y);
				z=ParseInt(sub[2],y);
					if (x<0||x>MX+1||y<0||y>MY+1||z<0||z>MZ+1) {
						cout <<"Requested output position is out of bounds." << endl;
						return 0;
					} else return X[JX*x+JY*y+z];
				} else  cout <<"Request for profile output does not contain the expected coordinate." << endl;
			return 0;
			break;
		default:
			break;
	}
	return 0;
}



bool Lattice::PrepareForCalculations(void) {
NAMICS_DBG("PrepareForCalculations in lattice" << endl);	bool success=true;
	return success;
}

void Lattice::push(string s, Real X) {
NAMICS_DBG("push (Real) in lattice " << endl); Reals.push_back(s);
	Reals_value.push_back(X);
}
void Lattice::push(string s, int X) {
NAMICS_DBG("push (int) in lattice " << endl); ints.push_back(s);
	ints_value.push_back(X);
}
void Lattice::push(string s, bool X) {
NAMICS_DBG("push (bool) in lattice " << endl); bools.push_back(s);
	bools_value.push_back(X);
}
void Lattice::push(string s, string X) {
NAMICS_DBG("push (string) in lattice " << endl); strings.push_back(s);
	strings_value.push_back(X);
}

Real* Lattice::GetPointer(string s,int &SIZE) {
NAMICS_DBG("GetPointer for lattice " + name << endl);	vector<string> sub;
	SIZE=M;
	In->split(s,';',sub);
	if (sub[0]=="profile" && sub[1]=="0") return L;
	if (sub[0]=="vector") {}
	return NULL;
}

int* Lattice::GetPointerInt(string s,int &SIZE) {
NAMICS_DBG("GetPointerInt for lattice " + name << endl);	vector<string> sub;
	SIZE=M;
	In->split(s,';',sub);
	if (sub[0]=="array"){//get with sub[1] the number and put the pointer to integer array in return.
	}

	return NULL;
}

void Lattice::PushOutput() {
NAMICS_DBG("PushOutput in lat " << endl); strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	string mirror="mirror";
	string periodic="periodic";
	push("geometry",geometry);
	push("gradients",gradients);
	if (offset_first_layer>0) push("offset_first_layer",offset_first_layer);
	push("volume",volume);
	push("accessible volume",Accesible_volume);
	string LatticeType;
	if (lattice_type == simple_cubic) LatticeType="simple_cubic"; else LatticeType="hexagonal";
	push("lattice_type",LatticeType);
	push("bond_length",bond_length);
	push("FJC_choices",FJC);
	string s="profile;0"; push("L",s);

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


int Lattice::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
NAMICS_DBG("GetValue (long)  in lattice " << endl);
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
