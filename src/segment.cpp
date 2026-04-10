#include "segment.h"
#include "io_utils.h"

Segment::Segment(const Input* In_,Lattice* Lat_, std::string name_) {
	In=In_; name=name_;
NAMICS_DBG("Segment constructor" + name << std::endl);
	lat=Lat_;
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
	u.assign(M * ns, 0);
	u_ext.assign(M, 0);
	phi.assign(M, 0);
	MASK.assign(M, 0);
	alpha.assign(M * ns, 0);
	phi_state.assign(M * ns, 0);
	G1.assign(M, 0);
	phi_side.assign(M * ns, 0);
	if (!ParseFreedoms()) std::cout <<"errors occurred.... progress uncertain...." << std::endl;

	all_segment=true;
}

bool Segment::ParseFreedoms() {
NAMICS_DBG("ParseFreedoms " << std::endl);
	const auto& parameters = In->Parameters("mon", name, start);
	MASK.assign(lat->M, 0);
	phibulk = 0;
	frozen_at_bound = -1;
	if (freedom == "free") return true;
	if (freedom == "pinned" && parameters.contains("frozen_range")) {
		std::cout << "For mon " << name << ", you should exclusively combine freedom : pinned with pinned_range" << std::endl;
		return false;
	}
	if (freedom == "frozen" && parameters.contains("pinned_range")) {
		std::cout << "For mon " << name << ", you should exclusively combine freedom : frozen with frozen_range" << std::endl;
		return false;
	}
	const char* key = freedom == "pinned" ? "pinned_range" : "frozen_range";
	const auto spec = parameters.find(key);
	if (spec == parameters.end()) {
		std::cout << "For mon " << name << ", you should provide '" << key << "'" << std::endl;
		return false;
	}
	const int nx = lat->MX / lat->fjc;
	const int ny = lat->gradients >= 2 ? lat->MY / lat->fjc : 1;
	const int nz = lat->gradients >= 3 ? lat->MZ / lat->fjc : 1;
	std::vector<std::vector<int>> coordinates;
	const auto& spec_value = spec.value();
	const auto values = spec_value.find("coordinates");
	if (!spec_value.is_object() || values == spec_value.end() || !values->is_array()) {
		std::cout << "mon " << name << " expects '" << key << "' as an object with a 'coordinates' array." << std::endl;
		return false;
	}
	try {
		coordinates = values->get<std::vector<std::vector<int>>>();
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in mon '" << name << "' for '" << key << "': " << error.what() << std::endl;
		return false;
	}
	for (const auto& point : coordinates) {
		if (static_cast<int>(point.size()) == lat->gradients) continue;
		std::cout << "mon " << name << " has an invalid coordinate in '" << key << "'." << std::endl;
		return false;
	}
	for (const auto& point : coordinates) {
		const int x = point[0];
		const int y = lat->gradients >= 2 ? point[1] : 0;
		const int z = lat->gradients >= 3 ? point[2] : 0;
		if (freedom == "pinned") {
			if (x < 1 || x > nx ||
			    (lat->gradients >= 2 && (y < 1 || y > ny)) ||
			    (lat->gradients >= 3 && (z < 1 || z > nz))) {
				std::cout << "For mon " << name << ", '" << key << "' contains a point outside the lattice." << std::endl;
				return false;
			}
		} else {
			const auto mark_surface = [&](int coordinate, int max_value, int low_bc, int high_bc) {
				if (coordinate == 0) {
					if (lat->BC[low_bc] != "surface") return false;
					frozen_at_bound = low_bc;
				} else if (coordinate == max_value + 1) {
					if (lat->BC[high_bc] != "surface") return false;
					frozen_at_bound = high_bc;
				} else if (coordinate < 0 || coordinate > max_value + 1) {
					return false;
				}
				return true;
			};
			if (!mark_surface(x, nx, 0, 3) ||
			    (lat->gradients >= 2 && !mark_surface(y, ny, 1, 4)) ||
			    (lat->gradients >= 3 && !mark_surface(z, nz, 2, 5))) {
				std::cout << "For mon " << name << ", '" << key << "' contains a point outside the accessible frozen range." << std::endl;
				return false;
			}
		}
		const int x0 = x * lat->fjc;
		const int x1 = (x + 1) * lat->fjc - 1;
		const int y0 = lat->gradients >= 2 ? y * lat->fjc : 0;
		const int y1 = lat->gradients >= 2 ? (y + 1) * lat->fjc - 1 : 0;
		const int z0 = lat->gradients >= 3 ? z * lat->fjc : 0;
		const int z1 = lat->gradients >= 3 ? (z + 1) * lat->fjc - 1 : 0;
		switch (lat->gradients) {
			case 1:
				for (int xx = x0; xx <= x1; ++xx) MASK[lat->P(xx)] = 1;
				break;
			case 2:
				for (int xx = x0; xx <= x1; ++xx) for (int yy = y0; yy <= y1; ++yy) MASK[lat->P(xx, yy)] = 1;
				break;
			case 3:
				for (int xx = x0; xx <= x1; ++xx) for (int yy = y0; yy <= y1; ++yy) for (int zz = z0; zz <= z1; ++zz) MASK[lat->P(xx, yy, zz)] = 1;
				break;
			default:
				return false;
		}
	}
	return true;
}

bool Segment::LoadExternalPotential() {
	std::fill(u_ext.begin(), u_ext.end(), 0);

	const auto& parameters = In->Parameters("mon", name, start);
	const std::string external_potential_filename = parameters.value("external_potential_filename", std::string{});
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

void Segment::PrepareForCalculations(std::span<const Real> KSAM, bool first_time) {
NAMICS_DBG("PrepareForCalcualtions in Segment " +name << std::endl);

	int M=lat->M;
	const auto& parameters = In->Parameters("mon", name, start);
	phibulk=0;
	if (freedom=="frozen") {
		std::copy_n(MASK.begin(), M, phi.begin());
	} else {
		std::fill(phi.begin(), phi.end(), 0);
	}

	if (!parameters.value("external_potential_filename", std::string{}).empty() && first_time) {
		if (!LoadExternalPotential()) return;
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
}

bool Segment::CheckInput(int start_) {
NAMICS_DBG("CheckInput in Segment " + name << std::endl);
	bool success;
	start=start_;
	unique=true;
	seg_nr_of_copy=-1;
	state_nr_of_copy=-1;
	ns=1;
	std::vector<std::string>options;

	fixedPsi0=false;
	const auto& parameters = In->Parameters("mon", name, start);
	static const std::vector<std::string> keys = {
		"freedom", "valence", "epsilon", "e.psi0/kT",
		"pinned_range", "frozen_range",
		"external_potential_filename", "var_pos", "set_equal_to"
	};
	chi.clear();
	success = true;
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		if (it.key().rfind("chi_", 0) == 0) {
			const std::string target = it.key().substr(4);
			if (ContainsValue(In->MonList, target) || ContainsValue(In->StateList, target)) {
				try {
					const Real chi_value = it.value().get<Real>();
					if (target == name && chi_value != 0) {
						std::cout <<" chi value for chi("<<name<<","<<target<<") value ignored: set to zero!" << std::endl;
					}
				} catch (const nlohmann::json::exception& error) {
					success = false;
					std::cout << "Invalid json type in mon '" << name << "' for '" << it.key() << "': " << error.what() << std::endl;
				}
				continue;
			}
		}
		success = false;
		std::cout << "mon property '" << it.key() << "' is unknown. Use a standard mon key or chi_<mon/state name>." << std::endl;
	}
	if(success) {
		try {
		if (parameters.contains("var_pos")) var_pos=parameters.at("var_pos").get<int>();

		std::string copy_of;
		if (parameters.contains("set_equal_to")) {
			copy_of=parameters.at("set_equal_to").get<std::string>();
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
			freedom = parameters.value("freedom", std::string{"free"});
			if (!ContainsValue(options,freedom)) {
				std::cout << "Freedom: '"<< freedom  <<"' for mon " + name + " not recognized. "<< std::endl;
				std::cout << "Freedom choices: free, pinned, frozen " << std::endl; success=false;
			}

		if (freedom =="free") {
			if (parameters.contains("frozen_range") || parameters.contains("pinned_range")) {
				if (start==1) {
					success=false;
					std::cout <<"In mon " + name + " you should not combine 'freedom : free' with 'frozen_range' or 'pinned_range'." << std::endl;
				}
			}
		}

		valence =0;
		if (parameters.contains("valence")) {
			valence=parameters.at("valence").get<Real>();
			if (valence<-10 || valence > 10) std::cout <<"For mon " + name + " valence value out of range -10 .. 10. Default value used instead" << std::endl;
		}
		epsilon=80;
		if (parameters.contains("epsilon")) {
			if (copy_of.size()>0) std::cout <<"For segment " << name << "value for epsilon will be overwritten by the value of segment " << copy_of << std::endl;
			epsilon=parameters.at("epsilon").get<Real>();
			if (epsilon<1 || epsilon > 250) std::cout <<"For mon " + name + " relative epsilon value out of range 1 .. 250. Default value 80 used instead" << std::endl;
		}
		if (valence !=0) {
			if (lat->bond_length <1e-12 || lat->bond_length > 1e-8) {
				success=false;
				if (lat->bond_length==0) std::cout << "When there are charged segments, you should set the bond_length in lattice to a reasonable value, e.g. between 1e-102... 1e-8 m " << std::endl;
				else std::cout <<"Bond length is out of range: 1e-12..1e-8 m " << std::endl;
			}
		}
		if (parameters.contains("e.psi0/kT")) {
			PSI0=0;
			fixedPsi0=true;
			PSI0=parameters.at("e.psi0/kT").get<Real>();
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
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in mon '" << name << "': " << error.what() << std::endl;
			success = false;
		}
	}

	int length = state_name.size();
	if (length >0 && freedom == "frozen") {
		success=false;
		std::cout <<" When freedom = 'frozen' a 'mon' can not have multiple internal states; status violated for mon " << name << std::endl;
	}

	if (parameters.contains("external_potential_filename")) {
		if (parameters.at("external_potential_filename").get<std::string>()=="?") {
			success=false;
			std::cout <<"Provide a json file containing an 'external_potential' std::array for mon " << name << std::endl;
		}
	}

	MASK.assign(lat->M, 0);
	if (success) success=ParseFreedoms();
	MASK.clear();
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

void Segment::PushOutput() {
NAMICS_DBG("PushOutput for segment " + name << std::endl);
	int M = lat->M;
	const auto& parameters = In->Parameters("mon", name, start);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["freedom"] = freedom;
	OUTPUT["valence"] = valence;
	Real theta=0;
	theta = lat->WeightedSum(phi.data());
	OUTPUT["theta"] = theta;
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
							OUTPUT["Range_z"] = pz;
							[[fallthrough]];
					case 2 :
							py=(point%JX)/JY;
							OUTPUT["Range_y"] = py;
							[[fallthrough]];
					case 1 :
							px=point/JX;
							OUTPUT["Range_x"] = px;
				break;
				default :
				break;
			}
		}
	}

	if (freedom != "frozen" && freedom != "pinned") theta_exc=theta-lat->volume*phibulk; else theta_exc=theta;
	OUTPUT["theta_exc"] = theta_exc;
	OUTPUT["phibulk"] = phibulk;
	if (parameters.contains("external_potential_filename")) OUTPUT["external_potential_filename"] = parameters.at("external_potential_filename").get<std::string>();
	if (freedom != "frozen" && freedom != "pinned") {
		OUTPUT["var_pos"] = var_pos;
	}
	if (freedom=="free") {
		Real first_moment = 0;
		Real second_moment = 0;
		Real fluctuations = 0;
		if (theta_exc !=0) first_moment=lat->Moment(phi.data(),phibulk,1)/theta_exc;
		if (theta_exc !=0) second_moment=lat->Moment(phi.data(),phibulk,2)/theta_exc;
		if (second_moment !=0) RMS=std::pow(second_moment,0.5);
		OUTPUT["RMS"] = RMS;
		OUTPUT["1st_M_phi_z"] = first_moment;
		OUTPUT["2nd_M_phi_z"] = second_moment;
		fluctuations = (second_moment-first_moment*first_moment);
		if (fluctuations >0) fluctuations = std::sqrt(fluctuations); else fluctuations=0;
		OUTPUT["fluctuations"] = fluctuations;
	}
	if (ns>1) {
		state_theta.clear();
		for (int i=0; i<ns; i++){
			OUTPUT["alphabulk_" + state_name[i]] = state_alphabulk[i];
			OUTPUT["valence_" + state_name[i]] = state_valence[i];
			OUTPUT["phibulk_" + state_name[i]] = state_phibulk[i];
			theta=lat->WeightedSum(phi_state.data()+i*M);
			state_theta.push_back(theta);
			OUTPUT["theta_" + state_name[i]] = theta;
			OUTPUT["theta_exc_" + state_name[i]] = theta - lat->volume * state_phibulk[i];
		}
	}
	for (size_t i = 0; i < In->MonList.size() && i < chi.size(); i++) OUTPUT["chi_" + In->MonList[i]] = chi[i];
	for (size_t i = 0; i < In->StateList.size() && In->MonList.size() + i < chi.size(); i++) OUTPUT["chi_" + In->StateList[i]] = chi[In->MonList.size() + i];
	if (fixedPsi0) OUTPUT["Psi0"] = PSI0;
	if (freedom=="pinned" && parameters.contains("pinned_range")) OUTPUT["range"] = parameters.at("pinned_range");
	if (freedom=="frozen" && parameters.contains("frozen_range")) OUTPUT["range"] = parameters.at("frozen_range");
	OUTPUT["phi"] = {{"profile", 0}};
	OUTPUT["G1"] = {{"profile", 1}};
	int profile = 2;
	if (lat->gradients==3) {
		OUTPUT["phi[z]"] = {{"profile", profile++}};
	}
	if (ns==1) {
		OUTPUT["u"] = {{"profile", profile++}};
	}
	if (ns >1) {
		for (int i = 0; i < ns; i++) {
			OUTPUT["phi-" + state_name[i]] = {{"profile", profile++}};
		}
		for (int i = 0; i < ns; i++) {
			OUTPUT["alpha-" + state_name[i]] = {{"profile", profile++}};
		}
		for (int i = 0; i < ns; i++) {
			OUTPUT["u-" + state_name[i]] = {{"profile", profile++}};
		}
	}
}

std::span<Real> Segment::GetPointer(int profile) {
NAMICS_DBG("Get Pointer for segment " + name << std::endl);
	int M=lat->M;
	if (profile == 0) {
		if (freedom=="frozen") {
			std::copy_n(MASK.begin(), M, phi.begin());
		} else lat->set_bounds(phi.data());
		return phi;
	}
	if (profile == 1) return G1;
	int offset = 2;
	if (lat->gradients == 3 && profile == offset++) {
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
	if (ns==1 && profile == offset) return u;
	if (ns>1) {
		if (profile >= offset && profile < offset + ns) {
			return std::span<Real>(phi_state).subspan(static_cast<size_t>(profile - offset) * M, static_cast<size_t>(M));
		}
		offset += ns;
		if (profile >= offset && profile < offset + ns) {
			return std::span<Real>(alpha).subspan(static_cast<size_t>(profile - offset) * M, static_cast<size_t>(M));
		}
		offset += ns;
		if (profile >= offset && profile < offset + ns) {
			return std::span<Real>(u).subspan(static_cast<size_t>(profile - offset) * M, static_cast<size_t>(M));
		}
	}
	return {};
}
void Segment::UpdateValence(Real*g, std::span<Real> psi, std::span<Real> q, std::span<Real> eps,bool grad_epsilon) {
	if (fixedPsi0) {
		std::transform(psi.begin(), psi.end(), MASK.begin(), psi.begin(),
		               [this](Real p, int mask_value) { return (mask_value == 1) ? PSI0 : p; });
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



void Segment::PutAlpha(Real alpha) { //expected to replace other method with same name.
	Real sum_alpha=0;
	int n_s;
	if (ns==1) n_s=0; else n_s=ns;
	if (n_s==0) return;

	for (int i=0; i<n_s; i++) {
		if (ItState !=i && state_change[i]) state_alphabulk[i]=0;
	}


	state_alphabulk[ItState]*=alpha;

	for (int i=0; i<n_s; i++) {
		sum_alpha+=state_alphabulk[i];
	}
	for (int i=0; i<n_s; i++) {
		if (state_alphabulk[i]==0) state_alphabulk[i]=1.0-sum_alpha;
		if (state_alphabulk[i]<0 || state_alphabulk[i]>1) std::cout <<"In Segment::PutAlpha, alphabulk out of bounds " << std::endl;
	}
}
