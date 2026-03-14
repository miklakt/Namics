#include "segment.h"
#include "io_utils.h"

Segment::Segment(const Input* In_,Lattice* Lat_, std::string name_,int segnr,int N_seg) {
	In=In_; name=name_; n_seg=N_seg; seg_nr=segnr;
NAMICS_DBG("Segment constructor" + name << std::endl);
	lat=Lat_;
	KEYS.push_back("freedom");
	KEYS.push_back("valence");
	KEYS.push_back("epsilon");
	KEYS.push_back("e.psi0/kT");
	KEYS.push_back("pinned_range");
	KEYS.push_back("frozen_range");
	KEYS.push_back("pinned_filename");
	KEYS.push_back("frozen_filename");
	KEYS.push_back("external_potential_filename");
	KEYS.push_back("var_pos");
	KEYS.push_back("set_equal_to");
	var_pos=0;
	ns=1;
	all_segment=false;
	phibulk=0;
	freedom="free";
}
Segment::~Segment() {
NAMICS_DBG("Segment destructor " + name << std::endl);
	DeAllocateMemory();
}

void Segment::DeAllocateMemory(void){
NAMICS_DBG( "In Segment, Deallocating memory " + name << std::endl);
if (!all_segment) return;
	u.clear();
	u_ext.clear();
	phi.clear();
	MASK.clear();
	alpha.clear();
	ALPHA.clear();
	phi_state.clear();
	G1.clear();
	phi_side.clear();
	all_segment=false;
}

void Segment::AllocateMemory() {
NAMICS_DBG("Allocate Memory in Segment " + name << std::endl);
	DeAllocateMemory();
	int M=lat->M;
	ns=state_name.size(); if (ns==0) ns=1;
	r.fill(0);
	u.assign(M * ns, 0);
	u_ext.assign(M, 0);
	phi.assign(M, 0);
	MASK.assign(M, 0);
	alpha.assign(M * ns, 0);
	ALPHA.assign(M * ns, 0);
	phi_state.assign(M * ns, 0);
	G1.assign(M, 0);
	phi_side.assign(M * ns, 0);
	bool success=true;
	bool HMaskDone=false;
	success=ParseFreedoms(HMaskDone);


	if (freedom!="free"&& !HMaskDone) {
		r[0]*=lat->fjc; r[1]*=lat->fjc; r[2]*=lat->fjc;
		r[3]=(r[3]+1)*lat->fjc-1;r[4]=(r[4]+1)*lat->fjc-1;r[5]=(r[5]+1)*lat->fjc-1;
		lat->CreateMASK(MASK.data(),r.data(),P.data(),n_pos,block);
	}
	if (!success) std::cout <<"errors occurred.... progress uncertain...." << std::endl;

	all_segment=true;
}

bool Segment::ParseFreedoms(bool& HMaskDone) {
NAMICS_DBG("ParseFreedoms " << std::endl);
	bool success=true;
	if (freedom == "pinned") {
		std::string freedom_source;
		int fjc = lat->fjc;
		int n_layers_x=(lat->MX)/fjc;
		int n_layers_y=(lat->MY)/fjc;
		int n_layers_z=(lat->MZ)/fjc;


		if (GetValue("pinned_range").size()>0) {
			freedom_source="pinned_range";
			std::string p_range=GetValue("pinned_range");
			std::vector<std::string>sub;
			In->split(p_range,';',sub);
			int Lsub=sub.size();
			std::vector<std::string>xyz;
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				for (int kk=0; kk<Lxyz; kk++){
					if (xyz[kk]=="firstlayer") {
						if (kk==0) p_range ="firstlayer_x";
						if (kk==1) p_range ="firstlayer_y";
						if (kk==2) p_range ="firstlayer_z";
					}
					if (xyz[kk]=="lastlayer") {
						if (kk==0) p_range ="lastlayer_x";
						if (kk==1) p_range ="lastlayer_y";
						if (kk==2) p_range ="lastlayer_z";
					}
				}
			}
			if (p_range=="firstlayer_x") {
				if (lat->gradients==1) {p_range = "firstlayer;firstlayer"; }
				if (lat->gradients==2) {p_range = "firstlayer,1;firstlayer,";p_range.append(std::to_string(n_layers_y)); }
				if (lat->gradients==3) {p_range = "firstlayer,1,1;firstlayer,"; p_range.append(std::to_string(n_layers_y)).append(",").append(std::to_string(n_layers_z)); }
			}
			if (p_range=="lastlayer_x") {
				if (lat->gradients==1) {p_range = "lastlayer;lastlayer";}
				if (lat->gradients==2) {p_range = "lastlayer,1;lastlayer,"; p_range.append(std::to_string(n_layers_y)); }
				if (lat->gradients==3) {p_range = "lastlayer,1,1;lastlayer,"; p_range.append(std::to_string(n_layers_y)).append(",").append(std::to_string(n_layers_z)); }
			}
			if (p_range=="firstlayer_y") {
				if (lat->gradients==2) {p_range = "1,firstlayer;";p_range.append(std::to_string(n_layers_x)).append(",firstlayer"); }
				if (lat->gradients==3) {p_range = "1,firstlayer,1;";p_range.append(std::to_string(n_layers_x)).append(",firstlayer,").append(std::to_string(n_layers_z)); }
			}
			if (p_range=="lastlayer_y") {
				if (lat->gradients==2) {p_range = "1,lastlayer;";p_range.append(std::to_string(n_layers_x)).append(",lastlayer"); }
				if (lat->gradients==3) {p_range = "1,lastlayer,1;";p_range.append(std::to_string(n_layers_x)).append(",lastlayer,").append(std::to_string(n_layers_z)); }
			}
			if (p_range=="firstlayer_z") {
				if (lat->gradients==3) {p_range = "1,1,firstlayer;";p_range.append(std::to_string(n_layers_x)).append(",").append(std::to_string(n_layers_y)).append(",firstlayer"); }
			}
			if (p_range=="lastlayer_z") {
				if (lat->gradients==3) {p_range = "1,1,lastlayer;";p_range.append(std::to_string(n_layers_x)).append(",").append(std::to_string(n_layers_y)).append(",lastlayer"); }
			}
			phibulk=0;
			if (GetValue("frozen_range").size()>0 || GetValue("frozen_filename").size()>0) {
			std::cout<< "For mon :" + name + ", you should exclusively combine freedom : pinned with pinned_range or pinned_filename" << std::endl;  success=false;}
			if (GetValue("pinned_range").size()>0 && GetValue("pinned_filename").size()>0) {
				std::cout<< "For mon " + name + ", you can not combine pinned_range with 'pinned_filename' " <<std::endl; success=false;
			}
			if (GetValue("pinned_range").size()==0 && GetValue("pinned_filename").size()==0) {
				std::cout<< "For mon " + name + ", you should provide either pinned_range or pinned_filename " <<std::endl; success=false;
			}
			sub.clear();
			In->split(p_range,';',sub);
			p_range.clear();
			Lsub=sub.size();

			if (Lsub!=2) {
				std::cout <<"For mon " + name + ", the parsing of 'pinned_range' failed. Use x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, or x1;x2 for 3, 2, or 1  gradient computations, respectively. x, y and z can also be  keys: 'firstlayer', 'lastlayer'" << std::endl;
				success=false;
				return success;
			}

			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				if (Lxyz<1 || Lxyz>3){
					std::cout <<"For mon " + name + ", the parsing of 'pinned_range' failed. Number of coordinates should be 1, 2 or 3: e.g., x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, x1;x2 for 1, 2 or 3 gradients, respectively.  " << std::endl;
					success=false;
					return success;
				}
				for (int kk=0; kk<Lxyz; kk++) {
					if (xyz[kk]=="firstlayer") {
						p_range.append("1");
					} else if (xyz[kk]=="lastlayer") {
						if (kk==0) p_range.append(std::to_string(n_layers_x));
						if (kk==1) p_range.append(std::to_string(n_layers_y));
						if (kk==2) p_range.append(std::to_string(n_layers_z));
					} else if (xyz[kk]=="var_pos") {
						if (((kk==0) && (var_pos<1 || var_pos> n_layers_x)) || ((kk==1) && (var_pos<1 || var_pos> n_layers_y))  ||((kk==2) && (var_pos<1 || var_pos> n_layers_z))) {
							std::cout <<"In pinned_range, 'var_pos' entry is out of bounds..." << std::endl; success=false; return success;
						}
						p_range.append(std::to_string(var_pos));

					} else {
						int cor=ParseInt(xyz[kk],-1);
						if (((kk==0) && (cor <1 || cor > n_layers_x)) || ((kk==1) && (cor <1 || cor > n_layers_y))  ||((kk==2) && (cor <1 || cor > n_layers_z))) {
							std::cout <<" For mon " + name+ ", the 'pinned_range' is not parsed properly! Coordinates either out of bounds or keywords 'var_pos', 'firstlayer', 'lastlayer' were not found" << std::endl;
							success=false;
							return success;
						} else p_range.append(xyz[kk]);
					}
					if (kk<Lxyz-1) p_range.append(",");
				}
				if (k<Lsub-1) p_range.append(";");
			}

			P.clear();
			n_pos=0;
			if (success) success=lat->ReadRange(r.data(), P.data(), n_pos, block, p_range,var_pos,name,freedom_source);
			if (n_pos>0) {
				P.assign(n_pos, 0);
				if (success) success=lat->ReadRange(r.data(), P.data(), n_pos, block, p_range,var_pos,name,freedom_source);
			}
		}
		if (GetValue("pinned_filename").size()>0) { freedom_source="pinned";
			block=false;
			const std::string filename=GetValue("pinned_filename");
			P.clear();
			n_pos=0;
			if (success) success=lat->ReadRangeFile(filename,P.data(),n_pos,name,freedom_source);
			if (n_pos>0) {
				P.assign(n_pos, 0);
				if (success) success=lat->ReadRangeFile(filename,P.data(),n_pos,name,freedom_source);
			}
		}
	}

	if (freedom == "frozen") {
		std::string freedom_source;
		int n_layers_x=(lat->MX)/lat->fjc;
		int n_layers_y=(lat->MY)/lat->fjc;
		int n_layers_z=(lat->MZ)/lat->fjc;

		frozen_at_bound=-1;
		phibulk=0;
		if (GetValue("pinned_range").size()>0 || GetValue("pinned_filename").size()>0) {
		        std::cout<< "For mon " + name + ", you should exclusively combine 'freedom : frozen' with 'frozen_range' or 'frozen_filename'" << std::endl;  success=false;
		}
		if (GetValue("frozen_range").size()>0 && GetValue("frozen_filename").size()>0) {
			std::cout<< "For mon " + name + ", you can not combine 'frozen_range' with 'frozen_filename' " <<std::endl; success=false;
		}
		if (GetValue("frozen_range").size()==0 && GetValue("frozen_filename").size()==0) {
			std::cout<< "For mon " + name + ", you should provide either 'frozen_range' or 'frozen_filename' " <<std::endl; success=false;
		}
		if (GetValue("frozen_range").size()>0) {
			freedom_source="frozen_range";
			std::string f_range=GetValue("frozen_range");
			std::vector<std::string>sub;
			In->split(f_range,';',sub);
			int Lsub=sub.size();
			std::vector<std::string>xyz;
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				for (int kk=0; kk<Lxyz; kk++){
					if (xyz[kk]=="lowerbound") {
						if (kk==0) f_range ="lowerbound_x";
						if (kk==1) f_range ="lowerbound_y";
						if (kk==2) f_range ="lowerbound_z";
					}
					if (xyz[kk]=="upperbound") {
						if (kk==0) f_range ="upperbound_x";
						if (kk==1) f_range ="upperbound_y";
						if (kk==2) f_range ="upperbound_z";
					}
				}
			}

			if (f_range=="lowerbound_x") {
				if (lat->gradients==1) {f_range = "lowerbound;lowerbound"; frozen_at_bound=0;}
				if (lat->gradients==2) {f_range = "lowerbound,1;lowerbound,"; f_range.append(std::to_string(n_layers_y)); frozen_at_bound=0;}
				if (lat->gradients==3) {f_range = "lowerbound,1,1;lowerbound"; f_range.append(std::to_string(n_layers_y)).append(",").append(std::to_string(n_layers_z)); frozen_at_bound=0;}
			}
			if (f_range=="upperbound_x") {
				if (lat->gradients==1) {f_range = "upperbound;upperbound";frozen_at_bound=3;}
				if (lat->gradients==2) {f_range = "upperbound,1;upperbound,"; f_range.append(std::to_string(n_layers_y)); frozen_at_bound=3;}
				if (lat->gradients==3) {f_range = "upperbound,1,1;upperbound"; f_range.append(std::to_string(n_layers_y)).append(",").append(std::to_string(n_layers_z)); frozen_at_bound=3;}
			}
			if (f_range=="lowerbound_y") {
				if (lat->gradients==2) {f_range = "1,lowerbound;";f_range.append(std::to_string(n_layers_x)).append(",lowerbound"); frozen_at_bound=1;}
				if (lat->gradients==3) {f_range = "1,lowerbound,1;";f_range.append(std::to_string(n_layers_x)).append(",lowerbound,").append(std::to_string(n_layers_z)); frozen_at_bound=1;}
			}
			if (f_range=="upperbound_y") {
				if (lat->gradients==2) {f_range = "1,upperbound;";f_range.append(std::to_string(n_layers_x)).append(",upperbound");  frozen_at_bound=4;}
				if (lat->gradients==3) {f_range = "1,upperbound,1;";f_range.append(std::to_string(n_layers_x)).append(",upperbound,").append(std::to_string(n_layers_z)); frozen_at_bound=4;}
			}

			if (f_range=="lowerbound_z") {
				if (lat->gradients==3) {f_range = "1,1,lowerbound;";f_range.append(std::to_string(n_layers_x)).append(",").append(std::to_string(n_layers_y)).append(",lowerbound"); frozen_at_bound=2; }
			}
			if (f_range=="upperbound_z") {
				if (lat->gradients==3) {f_range = "1,1,upperbound;";f_range.append(std::to_string(n_layers_x)).append(",").append(std::to_string(n_layers_y)).append(",upperbound"); frozen_at_bound=5;}
			}

			sub.clear();
			In->split(f_range,';',sub);
			f_range.clear();

			Lsub=sub.size();
			if (Lsub!=2) {
				std::cout <<"For mon " + name + ", the parsing of 'frozen_range' failed. Use x1,y1,z1;x2,y2,z2 in 3 gradients, x1,y1;x2,y2 in two gradients x1;x2 for one gradient computations. x, y and z can also be  keys: 'firstlayer', 'lastlayer', 'lowerbound', or 'upperbound'." << std::endl;
				std::cout <<"Alternatively - you can try a single keyword such as 'lowerbound', 'lowerbound_x', 'lowerbound_y', 'lowerbound_z', 'upperbound', 'upperbound_x', 'upperbound_y', 'upperbound_z'" << std::endl;
				success=false;
				return success;
			}
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				if (Lxyz<1 || Lxyz>3){
					std::cout <<"For mon " + name + ", the parsing of 'frozen_range' failed. Number of coordinates should be 1, 2 or 3: e.g., x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, x1;x2 for 1, 2 or 3 gradients, respectively.  " << std::endl;
					success=false;
					return success;
				}

				for (int kk=0; kk<Lxyz; kk++) {
					if (xyz[kk]=="firstlayer")
						f_range.append("1");
					if (xyz[kk]=="lastlayer") {
						if (kk==0) f_range.append(std::to_string(n_layers_x));
						if (kk==1) f_range.append(std::to_string(n_layers_y));
						if (kk==2) f_range.append(std::to_string(n_layers_z));
					}

					if (xyz[kk]=="lowerbound") {
						if ( (kk==0 && lat->BC[0] !="surface") ||(kk==1 && lat->BC[1] !="surface") ||(kk==2 && lat->BC[2] !="surface") ) {
							std::cout<<"In lattice you need boundary condition 'surface' in combination with frozen_range containing 'lowerbound' " << std::endl; success=false;
							return success;
						}
						f_range.append("0");
					}
					if (xyz[kk]=="upperbound") {
						if ( (kk==0 && lat->BC[3] !="surface") ||(kk==1 && lat->BC[4] !="surface") ||(kk==2 && lat->BC[5] !="surface") ) {
							std::cout<<"In lattice you need boundary condition 'surface' in combination with frozen_range containing 'upperbound' " << std::endl; success=false;
							return success;
						}
						if (kk==0) f_range.append(std::to_string(n_layers_x+1));
						if (kk==1) f_range.append(std::to_string(n_layers_y+1));
						if (kk==2) f_range.append(std::to_string(n_layers_z+1));
					}
					if (xyz[kk]=="var_pos") {
						if (((kk==0) && (var_pos<0 || var_pos> n_layers_x+1)) || ((kk==1) && (var_pos<0 || var_pos> n_layers_y+1))  ||((kk==2) && (var_pos<0 || var_pos> n_layers_z+1))) {
							std::cout <<"In frozen_range, 'var_pos' entry is out of bounds..." << std::endl; success=false; return success;
						}

						f_range.append(std::to_string(var_pos)); //check if var_pos is within lattice-range not implemented.
					}

					if (xyz[kk]!="firstlayer" && xyz[kk]!="lastlayer" && xyz[kk]!="lowerbound" && xyz[kk]!="upperbound" && xyz[kk]!="var_pos") {
						int cor=ParseInt(xyz[kk],-1);
						if ((kk==0 && (cor <0 || cor > n_layers_x+1)) || (kk==1 && (cor <0 || cor > n_layers_y+1))  ||(kk==2 && (cor <0 || cor > n_layers_z+1))) {
							std::cout <<" For mon " + name+ ", the 'frozen_range' is not parsed properly! Coordinates either out of bounds or keywords  'var_pos', 'firstlayer', 'lastlayer', 'lowerbound', 'upperbound' were not found" << std::endl;
							success=false;
							return success;
						} else {
							if (kk==0 && cor ==0 ) {
								if (lat->BC[0] !="surface") {std::cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << std::endl; success=false; return success;}
								frozen_at_bound=0;
							}
							if (kk==0 && cor == n_layers_x+1 ) {
								if (lat->BC[3] !="surface") {std::cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << std::endl; success=false; return success;}
								frozen_at_bound=3;
							}
							if (kk==1 && cor ==0 ) {
								if (lat->BC[1] !="surface") {std::cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << std::endl; success=false; return success;}
								frozen_at_bound=1;
							}
							if (kk==1 && cor == n_layers_y+1 ) {
								if (lat->BC[4] !="surface") {std::cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << std::endl; success=false; return success;}
								frozen_at_bound=4;
							}
							if (kk==2 && cor ==0 ) {
								if (lat->BC[2] !="surface") {std::cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << std::endl; success=false; return success;}
								frozen_at_bound=2;
							}
							if (kk==2 && cor == n_layers_z+1 ) {
								if (lat->BC[5] !="surface") {std::cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << std::endl; success=false; return success;}
								frozen_at_bound=5;
							}
							f_range.append(xyz[kk]);
						}
					}
					if (kk<Lxyz-1) f_range.append(",");
				}
				if (k<Lsub-1) f_range.append(";");
			}

			P.clear();
			n_pos=0;
			success=lat->ReadRange(r.data(), P.data(), n_pos, block, f_range,var_pos,name,freedom_source);
			if (n_pos>0) {
				P.assign(n_pos, 0);
				success=lat->ReadRange(r.data(), P.data(), n_pos, block, f_range,var_pos,name,freedom_source);
			}
		}
		if (GetValue("frozen_filename").size()>0) { freedom_source="frozen";
			block=false;
			const std::string filename=GetValue("frozen_filename");
			P.clear();
			n_pos=0;
			if (success) success=lat->ReadRangeFile(filename,P.data(),n_pos,name,freedom_source);
			if (n_pos>0) {
				P.assign(n_pos, 0);
				if (success) success=lat->ReadRangeFile(filename,P.data(),n_pos,name,freedom_source);
			}
		}
	}

	return success;
}

Real Segment::PinnedVolume() {
	int M=lat->M;
	Real volume=0;
	Real VOLUME=0;
	if (freedom !="pinned") return volume;
	if (lat->geometry=="planar") {
		(VOLUME) = 0; for (int __i = 0; __i < (M); ++__i) (VOLUME) += (MASK)[__i]; volume=1.0*VOLUME;
	} else {
		for (int i=0;i<M; i++) volume += MASK[i]*lat->L[i];
	}
	return volume/lat->fjc;
}

bool Segment::LoadExternalPotential() {
	std::fill(u_ext.begin(), u_ext.end(), 0);

	const std::string external_potential_filename = GetValue("external_potential_filename");
	if (external_potential_filename.size()==0) return true;
	const std::string resolved_external_potential_filename = In->ResolvePath(external_potential_filename);
	std::vector<Real> external_potential;
	if (!io::ReadExternalPotentialJson(resolved_external_potential_filename, external_potential)) {
		return false;
	}
	int expected = lat->MX;
	if (lat->gradients == 2) expected = lat->MX * lat->MY;
	if (lat->gradients == 3) expected = lat->MX * lat->MY * lat->MZ;
	if (static_cast<int>(external_potential.size()) != expected) {
		std::cout << "Inputfile " << resolved_external_potential_filename << " has " << external_potential.size()
		     << " values for 'external_potential', expected " << expected << " for mon " << name << std::endl;
		return false;
	}
	// Flattened json profile order matches Output::WriteOutput:
	// x-major in 1D, x/y-major in 2D, x/y/z-major in 3D.
	int pos = 0;
	if (lat->gradients == 1) {
		for (int x=1; x<=lat->MX; x++) u_ext[x] = external_potential[pos++];
	} else if (lat->gradients == 2) {
		for (int x=1; x<=lat->MX; x++) for (int y=1; y<=lat->MY; y++) u_ext[lat->P(x,y)] = external_potential[pos++];
	} else if (lat->gradients == 3) {
		for (int x=1; x<=lat->MX; x++) for (int y=1; y<=lat->MY; y++) for (int z=1; z<=lat->MZ; z++) u_ext[lat->P(x,y,z)] = external_potential[pos++];
	} else {
		std::cout << "Unsupported number of gradients for external_potential_filename in mon " << name << std::endl;
		return false;
	}
	return true;
}

bool Segment::PrepareForCalculations(std::span<const Real> KSAM, bool first_time) {
NAMICS_DBG("PrepareForCalcualtions in Segment " +name << std::endl);

	int M=lat->M;

	bool success=true;
	phibulk=0;
	if (freedom=="frozen") {
		std::copy_n(MASK.begin(), M, phi.begin());
	} else {
		std::fill(phi.begin(), phi.end(), 0);
	}

	if (GetValue("external_potential_filename").size()>0 && first_time) {
		success=LoadExternalPotential();
		if (!success) return false;
		if (ns==1) {
			for (int __i = 0; __i < (M); ++__i) (u)[__i] += (u_ext)[__i];
		} else {
			for (int i=0; i<ns; i++) {
				for (int __i = 0; __i < (M); ++__i) (u.data()+i*M)[__i] += (u_ext)[__i];
			}
		}
	}

	if (ns==1) {
		lat->set_bounds(u.data());
		for (int __i = 0; __i < (M); ++__i) (G1)[__i] = std::exp(-(u)[__i]);
	} else {
		std::fill(G1.begin(), G1.end(), 0);
		for (int i=0; i<ns; i++) {
			lat->set_bounds(u.data()+M*i);
			for (int __i = 0; __i < (M); ++__i) (alpha.data()+M*i)[__i] = std::exp(-(u.data()+M*i)[__i]);
			for (int __i = 0; __i < (M); ++__i) (alpha.data()+M*i)[__i] *= (state_alphabulk[i]);
			for (int __i = 0; __i < (M); ++__i) (G1)[__i] += (alpha.data()+M*i)[__i];
		}
		for (int i=0; i<ns; i++) for (int __i = 0; __i < (M); ++__i) (alpha.data()+i*M)[__i] = ((G1)[__i] != 0) ? ((alpha.data()+i*M)[__i] / (G1)[__i]) : 0;
	}

	if (freedom=="pinned") for (int __i = 0; __i < (M); ++__i) (G1)[__i] = (G1)[__i] * (MASK)[__i];
	if (freedom != "frozen") for (int __i = 0; __i < (M); ++__i) (G1)[__i] = (G1)[__i] * KSAM[__i];
	return success;
}

bool Segment::CheckInput(int start_) {
NAMICS_DBG("CheckInput in Segment " + name << std::endl);
	bool success;
	start=start_;
	block=false;
	unique=true;
	seg_nr_of_copy=-1;
	state_nr_of_copy=-1;
	ns=1;
	std::vector<std::string>options;
	n_pos=0;

	fixedPsi0=false;
	success = In->CheckParameters("mon",name,start, KEYS, PARAMETERS);
	if(success) {
		if (GetValue("var_pos").size()>0) var_pos=ParseInt(GetValue("var_pos"),0);

		std::string copy_of;
		if (GetValue("set_equal_to").size()>0) {
			copy_of=GetValue("set_equal_to");
			if (copy_of=="?") { success=false;
				std::cout <<" The following is expected: 'set_equal_to : segname' where 'segname' as a valid name of a segment. " << std::endl;
				std::cout <<" Use with caution. The 'epsilon' and all 'chi'-parameters of the segment " << name << " will be copied from the (indicated) segment " << std::endl;
				std::cout <<" Note that the value of 'valence' is not copied. You can/should still set the 'valence' of segment '" << name << "' uniquely." <<std::endl;
				std::cout <<" This has to do with the fact that 'electrostatics' is not effecting whether or not iteration variables are used for a particular segment. " <<std::endl;
			}
		}


			options.push_back("free");
			options.push_back("pinned");
			options.push_back("frozen");
			freedom="free";
			freedom = ParseString(GetValue("freedom"),"free");
			if (!In->InSet(options,freedom)) {
				std::cout << "Freedom: '"<< freedom  <<"' for mon " + name + " not recognized. "<< std::endl;
				std::cout << "Freedom choices: free, pinned, frozen " << std::endl; success=false;
			}

		if (freedom =="free") {
			if (GetValue("frozen_range").size()>0||GetValue("pinned_range").size()>0 ||
			GetValue("frozen_filename").size()>0 || GetValue("pinned_filename").size()>0) {
					if (start==1) {success=false; std::cout <<"In mon " + name + " you should not combine 'freedom : free' with 'frozen_range' or 'pinned_range' or corresponding filenames." << std::endl;
				}
			}
		}

		valence =0;
		if (GetValue("valence").size()>0) {
			valence=ParseReal(GetValue("valence"),0);
			if (valence<-10 || valence > 10) std::cout <<"For mon " + name + " valence value out of range -10 .. 10. Default value used instead" << std::endl;
		}
		epsilon=80;
		if (GetValue("epsilon").size()>0) {
			if (copy_of.size()>0) std::cout <<"For segment " << name << "value for epsilon will be overwritten by the value of segment " << copy_of << std::endl;
			epsilon=ParseReal(GetValue("epsilon"),80);
			if (epsilon<1 || epsilon > 250) std::cout <<"For mon " + name + " relative epsilon value out of range 1 .. 250. Default value 80 used instead" << std::endl;
		}
		if (valence !=0) {
			if (lat->bond_length <1e-12 || lat->bond_length > 1e-8) {
				success=false;
				if (lat->bond_length==0) std::cout << "When there are charged segments, you should set the bond_length in lattice to a reasonable value, e.g. between 1e-102... 1e-8 m " << std::endl;
				else std::cout <<"Bond length is out of range: 1e-12..1e-8 m " << std::endl;
			}
		}
		if (GetValue("e.psi0/kT").size()>0) {
			PSI0=0;
			fixedPsi0=true;
			PSI0=ParseReal(GetValue("e.psi0/kT"),0);
			if (PSI0!=0 && valence !=0) {
				success=false;
				std::cout <<"You can set only 'valence' or 'e.psi0/kT', but not both " << std::endl;
			}
			if (PSI0!=0 && freedom!="frozen") {
				success=false;
				std::cout <<"You can not set potential on segment that has not freedom 'frozen' " << std::endl;
			}
			if (PSI0 <-25 || PSI0 > 25) {
				success=false;
				std::cout <<"Value for dimensionless surface potentials 'e.psi0/kT' is out of range -25 .. 25. Recall the value of 1 at room temperature is equivalent to approximately 25 mV " << std::endl;
			}
		}
	}

	int length = state_name.size();
	if (length >0 && freedom == "frozen") {
		success=false;
		std::cout <<" When freedom = 'frozen' a 'mon' can not have multiple internal states; status violated for mon " << name << std::endl;
	}

	length=chi_name.size();

	Real Chi;
	for (int i=0; i<length; i++) {
		Chi=-999;
		const std::string chi_value = GetValue("chi_"+chi_name[i]);
		if (chi_value.size()>0) {
			Chi=ParseReal(chi_value,Chi);
			if (Chi==-999) {success=false; std::cout <<" chi value: chi("<<name<<","<<chi_name[i]<<") = "<<chi_value << "not valid." << std::endl; }
			if (name==chi_name[i] && Chi!=0) {if (Chi!=-999) std::cout <<" chi value for chi("<<name<<","<<chi_name[i]<<") = "<<chi_value << "value ignored: set to zero!" << std::endl; Chi=0;}

		}
		chi[i]=Chi;
	}

	if (GetValue("external_potential_filename").size()>0) {
		if (GetValue("external_potential_filename")=="?") {
			success=false;
			std::cout <<"Provide a json file containing an 'external_potential' std::array for mon " << name << std::endl;
		}
	}

	bool HMD=false;
	MASK.assign(lat->M, 0);
	P.clear();
	r.fill(0);
	if (success) success=ParseFreedoms(HMD);
	MASK.clear();
	P.clear();
	return success;
}

void Segment::SetPhiSide(){
NAMICS_DBG("SetPhiSide in Segment " + name << std::endl);
	int M=lat->M;
	if (ns==1) {
		lat->Side(phi_side.data(),phi.data(),M);
	} else {
		for (int i=0; i<ns; i++) {
			for (int __i = 0; __i < (M); ++__i) (phi_state.data()+i*M)[__i] = (alpha.data()+i*M)[__i] * (phi)[__i];
			lat->Side(phi_side.data()+i*M,phi_state.data()+i*M,M);
			state_phibulk[i]=phibulk*state_alphabulk[i];
		}
	}
}

void Segment::PutChiKEY(std::string new_name) {
NAMICS_DBG("PutChiKey " + name << std::endl);
	KEYS.push_back("chi_" + new_name);
	chi_name.push_back(new_name);
	chi.push_back(-999);
}

std::string Segment::GetValue(std::string parameter) {
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

void Segment::push(std::string s, Real X) {
NAMICS_DBG("Push in Segment (Real) " + name << std::endl);
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Segment::push(std::string s, int X) {
NAMICS_DBG("Push in Segment (int) " + name << std::endl);
	ints.push_back(s);
	ints_value.push_back(X);
}
void Segment::push(std::string s, bool X) {
NAMICS_DBG("Push in Segment (bool) " + name << std::endl);
	bools.push_back(s);
	bools_value.push_back(X);
}
void Segment::push(std::string s, std::string X) {
NAMICS_DBG("Push in Segment (std::string) " + name << std::endl);
	strings.push_back(s);
	strings_value.push_back(X);
}
void Segment::PushOutput() {
NAMICS_DBG("PushOutput for segment " + name << std::endl);
	int M = lat->M;

	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	ints.clear();
	ints_value.clear();
	Reals.clear();
	Reals_value.clear();
	push("freedom",freedom);
	push("valence",valence);
	Real theta=0;
	theta = lat->WeightedSum(phi.data());
	push("theta",theta);
	Real theta_exc=0;
	Real RMS=0;
	if (freedom == "frozen" || freedom == "pinned") {
		Real num_of_points;
		(num_of_points) = 0; for (int __i = 0; __i < (M); ++__i) (num_of_points) += (MASK)[__i];
		if (num_of_points==1) {
			int px=0,py=0,pz=0;
			int gradients=lat->gradients;
			int point=0;
			int JX=lat->JX;
			int JY=lat->JY;
			for (int i=0; i<M; i++) if (MASK[i]==1) point =i;
				switch (gradients)  {
					case 3 :
							pz=(point%JX)%JY;
							push("Range_z",pz);
							[[fallthrough]];
					case 2 :
							py=(point%JX)/JY;
							push("Range_y",py);
							[[fallthrough]];
					case 1 :
							px=point/JX;
							push("Range_x",px);
				break;
				default :
				break;
			}
		}
	}

	if (freedom != "frozen" && freedom != "pinned") theta_exc=theta-lat->volume*phibulk; else theta_exc=theta;
	push("theta_exc",theta_exc);
	push("phibulk",phibulk);
	if (GetValue("external_potential_filename").size()>0) push("external_potential_filename",GetValue("external_potential_filename"));
	if (freedom != "frozen" && freedom != "pinned") {
		push("var_pos",var_pos);
	}
	if (freedom=="free") {
		Real first_moment = 0;
		Real second_moment = 0;
		Real fluctuations = 0;
		if (theta_exc !=0) first_moment=lat->Moment(phi.data(),phibulk,1)/theta_exc;
		if (theta_exc !=0) second_moment=lat->Moment(phi.data(),phibulk,2)/theta_exc;
		if (second_moment !=0) RMS=std::pow(second_moment,0.5);
		push("RMS",RMS);
		push("1st_M_phi_z",first_moment);
		push("2nd_M_phi_z",second_moment);
		fluctuations = (second_moment-first_moment*first_moment);
		if (fluctuations >0) fluctuations = std::sqrt(fluctuations); else fluctuations=0;
		push("fluctuations",fluctuations);
	}
	if (ns>1) {
		state_theta.clear();
		for (int i=0; i<ns; i++){
			push("alphabulk_"+state_name[i],state_alphabulk[i]);
			push("valence_"+state_name[i],state_valence[i]);
			push("phibulk_"+state_name[i],state_phibulk[i]);
			theta=lat->WeightedSum(phi_state.data()+i*M);
			state_theta.push_back(theta);
			push("theta_"+state_name[i],theta);
			push("theta_exc_"+state_name[i],theta-lat->volume*state_phibulk[i]);
		}
	}
	int length=chi_name.size();
	for (int i=0; i<length; i++) push("chi_"+chi_name[i],chi[i]);
	if (fixedPsi0) push("Psi0",PSI0);
	if (freedom=="pinned") push("range",GetValue("pinned_range"));
	if (freedom=="frozen") push("range",GetValue("frozen_range"));
	std::string profile="profile;0"; push("phi",profile);

	profile="profile;1"; push("G1",profile);
	if (lat->gradients==3) {
		profile="profile;2"; push("phi[z]",profile);

	}
	int k=1;
	std::string s;
	std::string str;
	if (ns >1) {
		for (int i=0; i<ns; i++) { k++;
			std::stringstream ss; ss<<k; str=ss.str();
			s="profile;"+str; push("phi-"+state_name[i],s);
		}
		for (int i=0; i<ns; i++) { k++;
			std::stringstream ss; ss<<k; str=ss.str();
			s="profile;"+str; push("alpha-"+state_name[i],s);
		}
		for (int i=0; i<ns; i++) { k++;
			std::stringstream ss; ss<<k; str=ss.str();
			s="profile;"+str; push("u-"+state_name[i],s);
		}
	}


}

std::span<Real> Segment::GetPointer(std::string s) {
NAMICS_DBG("Get Pointer for segment " + name << std::endl);
	std::vector<std::string> sub;
	int M=lat->M;
	In->split(s,';',sub);
	if (sub[0]=="profile") {

	if (sub[1]=="0") {
		if (freedom=="frozen") {
			std::copy_n(MASK.begin(), M, phi.begin());
		} else lat->set_bounds(phi.data());
		return phi;
	}
	if (sub[1]=="1") return G1;
        if (sub[1]=="2") {
		int MX=lat->MX;
		int MY=lat->MY;
		int MZ=lat->MZ;
		int JX=lat->JX;
		int JY=lat->JY;
		Real Sum;
		std::fill(phi_side.begin(), phi_side.begin() + M, 0); //phi_side is reused because this std::array is no longer needed (hopefully....).
		for (int z=0; z<MZ; z++) {
			Sum=0;
			for (int x=1; x<MX+1; x++) for (int y=1; y<MY+1; y++)
				Sum +=phi[x*JX+y*JY+z];
			Sum /=MX*MY;
			phi_side[JX+JY+z]=Sum;
		}
		return phi_side;
	}
	if (ns>1) {
		for (int i=0; i<ns; i++) {
			std::stringstream ss; ss<<i+2; std::string str=ss.str();
			if (sub[1]==str) return std::span<Real>(phi_state).subspan(static_cast<size_t>(i * M), static_cast<size_t>(M));
		}
		for (int i=0; i<ns; i++) {
			std::stringstream ss; ss<<i+ns+2; std::string str=ss.str();
			if (sub[1]==str) return std::span<Real>(alpha).subspan(static_cast<size_t>(i * M), static_cast<size_t>(M));
		}
		for (int i=0; i<ns; i++) {
			std::stringstream ss; ss<<i+2*ns+2; std::string str=ss.str();
			if (sub[1]==str) return std::span<Real>(u).subspan(static_cast<size_t>(i * M), static_cast<size_t>(M));
		}
	}


	} else {//sub[0]=="std::vector" ..do not forget to set SIZE before returning the pointer.
	}
	return {};
}
std::span<int> Segment::GetPointerInt(std::string s) {
NAMICS_DBG("GetPointerInt for segment " + name << std::endl);
	std::vector<std::string> sub;
	In->split(s,';',sub);
	if (sub[0]=="std::array") {// set SIZE and return int pointer.
	}
	return {};
}

int Segment::GetValue(std::string prop,int &int_result,Real &Real_result,std::string &string_result){
NAMICS_DBG("GetValue long for segment " + name << std::endl);
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true"; else string_result="false";
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}
		i++;
	}
	return 0;
}
void Segment::UpdateValence(Real*g, std::span<Real> psi, std::span<Real> q, std::span<Real> eps,bool grad_epsilon) {
	int M=lat->M;
	if (fixedPsi0) {

		OverwriteC(psi.data(),MASK.data(),PSI0,M);
		lat->UpdateQ(g,psi.data(),q.data(),eps.data(),MASK.data(),grad_epsilon);
	}

}
int Segment::AddState(int id_,Real alphabulk,Real valence,bool fixed) {
NAMICS_DBG("AddState " << id_ <<" to seg " << name << std::endl);
	int length = state_name.size();
	int state_number=-1;
	bool found=false;
	int ID=id_;
	for (int k=0; k<length; k++) {
		if (state_id[k]==ID) {
			state_name[k]=In->StateList[ID];
			found=true;
			state_alphabulk[k]=alphabulk;
			state_valence[k]=valence;
			if (fixed) state_change[k]=false; else state_change[k]=true;
			state_number=k;
		}
	}
	if (!found) {
		state_id.push_back(ID);
		state_name.push_back(In->StateList[ID]);
		state_alphabulk.push_back(alphabulk);
		state_phibulk.push_back(0);
		state_valence.push_back(valence);
		if (fixed) state_change.push_back(false); else state_change.push_back(true);
		state_number=state_change.size()-1;
		length=In->StateList.size();
		for (int k=0; k<length; k++) {if (name==In->StateList[k]) state_nr.push_back(k);}
	}

	if (valence !=0) {
		if (lat->bond_length <1e-10 || lat->bond_length > 1e-8) {
			if (lat->bond_length==0) std::cout << "When there are charged states, you should set the bond_length in lattice to a reasonable value, e.g. between 1e-10 ... 1e-8 m " << std::endl;
			else std::cout <<"Bond length is out of range: 1e-10..1e-8 m " << std::endl;
		}
	}

	return state_number;
}



bool Segment::PutAlpha(Real alpha) { //expected to replace other method with same name.
	bool success=true;
	Real fixed_value=0;
	Real sum_alpha=0;
	int n_s;
	if (ns==1) n_s=0; else n_s=ns;
	if (n_s==0) return false;


	for (int i=0; i<n_s; i++) {
		if (!state_change[i]) fixed_value+=state_alphabulk[i];}
	for (int i=0; i<n_s; i++) {
		if (ItState !=i && state_change[i]) state_alphabulk[i]=0;
	}


	state_alphabulk[ItState]*=alpha;

	for (int i=0; i<n_s; i++) {
		sum_alpha+=state_alphabulk[i];
	}
	for (int i=0; i<n_s; i++) {
		if (state_alphabulk[i]==0) state_alphabulk[i]=1.0-sum_alpha;
		if (state_alphabulk[i]<0 || state_alphabulk[i]>1) {
			success=false; std::cout <<"In Segment::PutAlpha, alphabulk out of bounds " << std::endl;
		}
	}


	return success;
}

bool Segment::CanBeReached(int x0, int y0, int z0, int Ds) {
	int gradients=lat->gradients;
	int MX = lat->MX;
	int MY = lat->MY;
	int MZ = lat->MZ;
	int JX = lat->JX;
	int JY = lat->JY;
	int JZ = lat->JZ;
	int x=0,y=0,z=0;
	switch (gradients) {
		case 3:
			for (z=1; z<MZ+1; z++)
		case 2:
			for (y=1; y<MY+1; y++)
		case 1:
			for (x=1; x<MX+1; x++){
				if (MASK[x*JX+y*JY+z*JZ]==1) {
					if (std::abs(x-x0)+std::abs(y-y0)+std::abs(z-z0) < Ds) return true;
				}

			}
	}
	return false;
}
