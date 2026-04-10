#include "solve_scf.h"

Solve_scf::Solve_scf(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::span<const std::unique_ptr<State>> Sta_, std::span<const std::unique_ptr<Reaction>> Rea_, std::span<const std::unique_ptr<Molecule>> Mol_,System* Sys_,std::string name_) :
	name{name_}, In{In_}, Sys{Sys_}, Seg{Seg_}, lat{Lat_}, Mol{Mol_}, Sta{Sta_}, Rea{Rea_}
{
NAMICS_DBG("Constructor in Solve_scf " << std::endl);
	max_g = false; // compute g based on max error
	all=false;
	restart_DIIS =0;
}

Solve_scf::~Solve_scf() {
	DeAllocateMemory();
}

void Solve_scf :: DeAllocateMemory(){
NAMICS_DBG("DeAllocateMemory in Solve " << std::endl);
	if (!all) return;
	xx.clear();
	yy.clear();
	SIGN.clear();
	all=false;
NAMICS_DBG("exit for 'destructor' in Solve " << std::endl);

}

void Solve_scf::AllocateMemory() {
NAMICS_DBG("AllocateMemeory in Solve " << std::endl);
	if (all) DeAllocateMemory();
	int M=lat->M;
	iv = (Sys->ItMonList.size() + Sys->ItStateList.size())* M;
	if (Sys->charged) iv += M;
	xx.assign(iv, 0);
	int niv = In->ReactionList.size();
	if (niv>0) {
		yy.assign(niv, 0);
		SIGN.assign(niv, 1);
	}
	all=true;

	Sys->AllocateMemory();
}

bool Solve_scf::CheckInput(int start_) { start=start_;
NAMICS_DBG("CheckInput in Solve " << std::endl);
	pseudohessian =false;
	s_info=false;
	e_info=false;
	t_info=false;
	hs_info=false;
	i_info=1;
	hessian =false;
	bool success=true;
	control=proceed;
	solver=PSEUDOHESSIAN;
	SCF_method="pseudohessian";
	gradient=classical;
	residual=1;
	m=10;
	const auto& parameters = In->Parameters("newton", name, start);
	static const std::vector<std::string> keys = {
		"method", "e_info", "s_info", "i_info", "t_info", "hs_info",
		"iterationlimit", "tolerance", "stop_criterion", "deltamin", "deltamax",
		"linesearchlimit", "max_accuracy_for_hessian_scaling", "n_iterations_for_hessian",
		"small_alpha", "max_n_small_alpha", "min_accuracy_for_hessian", "m", "n_restart_DIIS"
	};
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		success = false;
		std::cout << "newton property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : keys) std::cout << item << std::endl;
	}
	if (success) {
		try {
		iterationlimit=parameters.value("iterationlimit",1000);
		if (iterationlimit < 0 || iterationlimit>1e6) {iterationlimit = 1000;}

		e_info=parameters.value("e_info",true);
		hs_info=parameters.value("hs_info",true);
		s_info=parameters.value("s_info",false);
		t_info=parameters.value("t_info",false);
		i_info=parameters.value("i_info",1);
		if (i_info == 0) {
		// We cannot divide by zero (see modulus statements in sfnewton), but this will probably be what the user means.
		std::cerr << "WARNING: i_info cannot be zero ! Defaulting to iterationlimit + 1."<< std::endl;
		i_info = iterationlimit+1;
		}
		deltamax=parameters.value("deltamax",0.1);
		if (deltamax < 0 || deltamax>100) {deltamax = 0.1;  std::cout << "Value of deltamax out of range 0..100, and value set to default value 0.1" <<std::endl; }
		deltamin=parameters.value("deltamin",0.0);
		if (deltamin < 0 || deltamin>100) {deltamin = deltamax/100000;  std::cout << "Value of deltamin out of range 0..100, and value set to default value deltamax/100000" <<std::endl; }
		tolerance=parameters.value("tolerance",1e-7);
		if (tolerance < 1e-16 ||tolerance>10) {tolerance = 1e-5;  std::cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<std::endl; }

			const std::string method_value = parameters.value("method", std::string{});
			if (method_value.empty()) {
				SCF_method="pseudohessian";
			} else if (method_value == "DIIS" || method_value == "pseudohessian" || method_value == "hessian" || method_value == "LBFGS") {
				SCF_method = method_value;
			} else {
				std::cout << "In 'solve_scf' the entry for 'method' not recognized: choose from:" << std::endl;
				std::cout << "DIIS" << std::endl;
				std::cout << "pseudohessian" << std::endl;
				std::cout << "hessian" << std::endl;
				std::cout << "LBFGS" << std::endl;
				success = false;
			}
		if (SCF_method=="hessian" || SCF_method=="pseudohessian") {
			if (SCF_method=="hessian") {pseudohessian=false; hessian=true; solver=HESSIAN;} else { pseudohessian=true; hessian=false; solver=PSEUDOHESSIAN;}
			samehessian=false;
			max_accuracy_for_hessian_scaling=parameters.value("max_accuracy_for_hessian_scaling",0.1);
			if (max_accuracy_for_hessian_scaling<1e-7 || max_accuracy_for_hessian_scaling>1) {
				std::cout <<"max_accuracy_for_hessian_scaling is out of range: 1e-7...1; default value 0.1 is used instead" << std::endl;
				max_accuracy_for_hessian_scaling=0.1;
			}
			minAccuracyForHessian=parameters.value("min_accuracy_for_hessian",0.5);
			if (minAccuracyForHessian<0 ||minAccuracyForHessian>1) {
				std::cout <<"min_accuracy_for_hessian is out of range: 0...0.1; default value 0 is used instead (no hessian computation)" << std::endl;
				minAccuracyForHessian=0;
			}
			n_iterations_for_hessian=parameters.value("n_iterations_for_hessian",iterationlimit+100);
			if (n_iterations_for_hessian<1 ) {
				std::cout <<" n_iterations_for_hessian setting must be larger than unity; hessian evaluations will not be done " << std::endl;
				n_iterations_for_hessian=iterationlimit+100;
			}
			maxNumSmallAlpha=parameters.value("max_n_small_alpha",50);
			if (maxNumSmallAlpha<10 ||maxNumSmallAlpha>1000) {
				std::cout <<" max_n_small_alpha is out of range: 10, ..., 100;  max_n_small_alpha is set to default: 50 " << std::endl;
				maxNumSmallAlpha=50;
			}
			deltamin=parameters.value("deltamin",0.0);
			if (deltamin <0 || deltamin>deltamax) {
				std::cout <<"deltamin is out of range; 0, ..., " << deltamax << "; deltamin value set to 0 " << std::endl;
				deltamin=0;
			}
			smallAlpha=parameters.value("small_alpha",0.00001);
			if (smallAlpha <0 || smallAlpha>1) {
				std::cout <<"small_alpha is out of range; 0, ..., 1; small_alpha value set to default: 1e-5 " << std::endl;
				smallAlpha=0.00001;
			}
		}
		if (SCF_method=="DIIS") {
			solver=diis;
			m=parameters.value("m",10);
			if (m < 0 ||m>100) {m=10;  std::cout << "Value of 'm' out of range 0..100, value set to default value 10" <<std::endl; }
			restart_DIIS=parameters.value("n_restart_DIIS",iterationlimit);
			if (restart_DIIS < 0 || restart_DIIS > iterationlimit*10) {
				restart_DIIS=iterationlimit; std::cout <<"Value of 'n_restart_DIIS' out of range 0 .. iterationlimit; value set to iterationlimit" << std::endl;
			}
			restart_DIIS -=restart_DIIS%m; std::cout <<"Restart DIIS set to " << restart_DIIS << std::endl;
		}
		if (SCF_method=="LBFGS") {
			solver=LBFGS;
			m=parameters.value("m",6);
			if (m < 0 ||m>1000) {m=6;  std::cout << "Value of 'm' out of range 0..1000, value set to default value 6" <<std::endl; }
		}
		const std::string stop_value = parameters.value("stop_criterion", std::string{});
		if (!stop_value.empty()) {
			if (stop_value != "norm_of_g" && stop_value != "max_of_element_of_|g|") {
				std::cout << "In newton the stop_criterion setting was not recognised" << std::endl;
				success = false;
			} else {
				stop_criterion = stop_value;
			}
			if(stop_value == "max_of_element_of_|g|") {
				max_g = true;
			}
		}
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in newton '" << name << "': " << error.what() << std::endl;
			success = false;
		}
	}
	return success;
}
void Solve_scf::PushOutput() {
NAMICS_DBG("PushOutput in  Solve " << std::endl);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["method"] = SCF_method;
	OUTPUT["m"] = m;
	OUTPUT["delta_max"] = deltamax;
	OUTPUT["residual"] = residual;
	OUTPUT["tolerance"] = tolerance;
	OUTPUT["iterations"] = iterations;
	OUTPUT["iterationlimit"] = iterationlimit;
	OUTPUT["stop_criterion"] = stop_criterion;
	if (pseudohessian || hessian) {
		OUTPUT["linesearchlimit"] = linesearchlimit;
		OUTPUT["max_accuracy_for_hessian_scaling"] = max_accuracy_for_hessian_scaling;
		OUTPUT["n_iteratons_for_hessian"] = n_iterations_for_hessian;
		OUTPUT["small_alpha"] = smallAlpha;
		OUTPUT["max_n_small_alpha"] = maxNumSmallAlpha;
		OUTPUT["min_accuracy_for_hessian"] = minAccuracyForHessian;
	}
	lat->PushOutput();
	int length = In->MonList.size();
	for (int i=0; i<length; i++) {
		Seg[i]->PushOutput();
	}
	length = In->MolList.size();
	for (int i=0; i<length; i++){
		Mol[i]->PushOutput();
	}
	length = In->StateList.size();
	for (int i=0; i<length; i++) Sta[i]->PushOutput();
	length = In->ReactionList.size();
	for (int i=0; i<length; i++) Rea[i]->PushOutput();
	Sys->PushOutput();
}

void Solve_scf::Copy(std::span<Real> x, std::span<const Real> X, int MX, int MY, int MZ, int fjc_old) {
	int mx=lat->MX;
	int my=lat->MY;
	int mz=lat->MZ;
	int jx=lat->JX;
	int jy=lat->JY;
	int i,j,k;
	int pos_i,pos_o;
	int fjc=lat->fjc;
	int JX=(MY+2*fjc_old)*(MZ+2*fjc_old);
	int JY=(MZ+2*fjc_old);


	switch (lat->gradients) {
		case 1:
			if (fjc==1 and fjc_old==1) {
				if (MY>0||MZ>0) {
					std::cout <<" Copy from more than one gradient to one gradient: (i) =(1,i) or (1,1,i) is used "<< std::endl;
				}
				if (MZ>0) { pos_i=JX+JY; pos_o=MZ+2;} else {if (MY>0) {pos_i=JX; pos_o=MY+2; } else { pos_i=0; pos_o=MX+2; } }
				for (i=0; i<mx+2*fjc; i++)  if (i<pos_o) x[i]=X[pos_i+i];
			} else {
				for (i=0; i<mx+2*fjc; i++) {
					x[i]=X[i];
				}
			}
			break;
		case 2:
			if (MY==0) {
				std::cout <<" Copy from one-gradient to two gradients: one-gradient (x,i)=(i) is used for all x " << std::endl;
				for (i=0; i<mx+2*fjc; i++)
				for (j=0; j<my+2*fjc; j++) if (j<MX+2*fjc_old) x[i*jx+j]=X[j];
			} else {
				if (MZ>0) {
					std::cout <<" Copy from three gradients to two gradients: (i,j)=(1,i,j) is used " <<std::endl;
					JX=(MY+2*fjc_old)*(MZ+2*fjc_old);
					JY=(MZ+2*fjc_old);
					for (i=0; i<mx+2*fjc; i++)
					for (j=0; j<my+2*fjc; j++) if (i<MY+2*fjc_old && j<MZ+2*fjc_old) x[i*jx+j]=X[JX+i*JY+j];
				} else {
					JX=(MY+2*fjc_old);
					for (i=0; i<mx+2*fjc; i++)
					for (j=0; j<my+2*fjc; j++) if (i<MX+2*fjc_old && j<MY+2*fjc_old) x[i*jx+j]=X[i*JX+j];
				}
			}
			break;
		case 3:
				if (MY==0) {
					std::cout <<"Copy from one gradient to three gradients: (x,y,i) = (i) is used for all x,y " << std::endl;
					for (i=0; i<mx+2*fjc; i++)
					for (j=0; j<my+2*fjc; j++)
					for (k=0; k<mz+2*fjc; k++) if (k<MX+2*fjc_old) x[i*jx+j*jy+k] = X[k];
				} else {
					if (MZ==0) {
						std::cout <<"Copy form two gradients to three: (x,i,j) = (i,j) for all x " << std::endl;
						JX=(MY+2*fjc_old);
						for (i=0; i<mx+2*fjc; i++)
						for (j=0; j<my+2*fjc; j++)
						for (k=0; k<mz+2; k++) if (j<MX+2*fjc_old && k<MY+2*fjc_old) x[i*jx+j*jy+k] = X[j*JX+k];
					} else {
						JX=(MY+2*fjc_old)*(MZ+2*fjc_old);
						JY=(MZ+2*fjc_old);
						for (i=0; i<mx+2*fjc; i++)
						for (j=0; j<my+2*fjc; j++)
						for (k=0; k<mz+2*fjc; k++) if (i<MX+2*fjc_old && j<MY+2*fjc_old && k<MZ+2*fjc_old) x[i*jx+j*jy+k] = X[i*JX+j*JY+k];
					}
				}
			break;
		default:
			break;
	}
}

void Solve_scf::Guess(std::span<const Real> X, std::vector<std::string> MONLIST, std::vector<std::string> STATELIST, bool CHARGED, int MX, int MY, int MZ,int fjc_old){
	NAMICS_DBG( "Guess in Solve" << std::endl);
	int M=lat->M;
	int m;
	if (MZ>0) {m=(MX+2)*(MY+2)*(MZ+2); } else { if (MY>0) { m=(MX+2*fjc_old)*(MY+2*fjc_old); } else {  m=(MX+2*fjc_old);}}

	int length_old_mon=MONLIST.size();
	int length_old_state=STATELIST.size();
	int length_new_mon=Sys->ItMonList.size();
	int length_new_state=Sys->ItStateList.size();
	auto xx_span = std::span<Real>(xx);
	for (int i = 0; i<length_old_mon; i++) {
		for (int j=0; j<length_new_mon; j++) {
			if (MONLIST[i]==Seg[Sys->ItMonList[j]]->name) {
				Copy(xx_span.subspan(static_cast<size_t>(M * j), static_cast<size_t>(M)),
				     X.subspan(static_cast<size_t>(i * m), static_cast<size_t>(m)),
				     MX,MY,MZ,fjc_old);
			}
		}
	}
	for (int i = 0; i<length_old_state; i++) {
		for (int j=0; j<length_new_state; j++) {
			if (STATELIST[i]==Sta[Sys->ItStateList[j]]->name) {
				Copy(xx_span.subspan(static_cast<size_t>(M * (j + length_new_mon)), static_cast<size_t>(M)),
				     X.subspan(static_cast<size_t>((i + length_old_mon) * m), static_cast<size_t>(m)),
				     MX,MY,MZ,fjc_old);
			}
		}
	}

	if (CHARGED && Sys->charged) {
		Copy(xx_span.subspan(static_cast<size_t>((length_new_mon + length_new_state) * M), static_cast<size_t>(M)),
		     X.subspan(static_cast<size_t>((length_old_mon + length_old_state) * m), static_cast<size_t>(m)),
		     MX,MY,MZ,fjc_old);
	}
}

class SCF_LBFGS
{
private:
    System* Sys;
    int iterations =0;
public:
    explicit SCF_LBFGS(System* Sys_) : Sys(Sys_) {}

	Real operator()(Vector& x_, Vector& g_)
    {
	int iv=x_.size();
	Sys->Classical_residual(std::span<const Real>(x_.data(), static_cast<size_t>(iv)),
	                        std::span<Real>(g_.data(), static_cast<size_t>(iv)),
	                        iterations);
	iterations++;
	return g_.norm();
	    }
};

bool Solve_scf::Solve(bool report_errors_) { //going SCF here
NAMICS_DBG("Solve in  Solve_scf " << std::endl);
	bool success = false;
	int niv = In->ReactionList.size();
	if (niv>0) {
		const iteration_method saved_solver = solver;
		const gradient_method saved_gradient = gradient;
		const inner_iteration_method saved_control = control;
		const bool saved_pseudohessian = pseudohessian;
		const bool saved_hessian = hessian;
		const bool saved_e_info = e_info;
		const bool saved_s_info = s_info;
		e_info = false;
		s_info = false;
		gradient = WEAK;
		control = super;
		pseudohessian = false;
		hessian = true;

		success = iterate(yy.data(), niv, 100, 1e-8, 1, 0.0000001, true);
		std::cout << iterations << " iterations to find alphabulk values. " <<std::endl;
		if (!success) std::cout <<"iteration for alphabulk values for internal states failed. Check eqns. " << std::endl;
		solver = saved_solver;
		gradient = saved_gradient;
		control = saved_control;
		pseudohessian = saved_pseudohessian;
		hessian = saved_hessian;
		e_info = saved_e_info;
		s_info = saved_s_info;
	}

	switch(solver) {
		case HESSIAN:
		case PSEUDOHESSIAN:
			success = iterate(xx.data(), iv, iterationlimit, tolerance, deltamax, deltamin, true);
		break;
		case diis:
			success = iterate_DIIS(xx.data(), iv, m, iterationlimit, tolerance, deltamax, restart_DIIS);
		break;
		case LBFGS:
			success = true;
			{
				SCF_LBFGS fun(Sys);
				LBFGSParam<Real> param;
				param.epsilon = tolerance;
				param.m = m;
				param.max_iterations = iterationlimit;
				if (deltamax > 0) param.max_step = deltamax;
				if (deltamin > 0 && deltamax > deltamin) param.min_step = deltamin;
				LBFGSSolver<Real> mysolver(param);
				Real fx = 0;
				std::cout <<std::endl <<"LBFGS has been notified" << std::endl;
				Vector x_vec = Eigen::Map<Vector>(xx.data(), iv);
				iterations = mysolver.minimize(fun, x_vec, fx);
				std::copy_n(x_vec.data(), iv, xx.begin());
				Real res = mysolver.final_grad_norm();
				std::cout <<std::endl <<"Problem solved: " << iterations << " iterations,  |g|: " << res <<  std::endl;
			}
		break;
	}
	if (success) Sys->FinalizeOutputs();
	return Sys->CheckResults(report_errors_);
}


void Solve_scf::residuals(Real* x, Real* g){
 NAMICS_DBG("residuals in Solve_scf " << std::endl);
	switch(gradient) {
		case WEAK:
			NAMICS_DBG("Residuals for weak iteration " << std::endl);
			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				if (Rea[i]->Sto.size()==3) Rea[i]->GuessAlpha();
			}
			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				if (Rea[i]->Sto.size()!=3) Rea[i]->GuessAlpha();
			}
			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				Rea[i]->PutAlpha(std::exp(x[i]));
			}


			std::fill_n(g, In->ReactionList.size(), 0);

			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				g[i]=SIGN[i]*Rea[i]->Residual_value();
			}
		break;
		default:
			NAMICS_DBG("Residuals in scf mode in Solve_scf " << std::endl);
			Sys->Classical_residual(std::span<const Real>(x, static_cast<size_t>(iv)),
			                        std::span<Real>(g, static_cast<size_t>(iv)),
			                        iterations);
		break;
	}
}

void Solve_scf::inneriteration(Real* g, Real* h, Real accuracy, Real& deltamax, Real ALPHA, int nvar) {
NAMICS_DBG("inneriteration in Solve_scf " << std::endl);
	residual=accuracy; // track the reported residual alongside the active solver accuracy.
	switch(control) {
		case super:
			for (int i=0; i<nvar; i++) { //this is to control the sing in the WEAK iteration.
				if (h[i+i*nvar]<0) {SIGN[i] = -1;  }
			}

		break;
		default:
			if (iterations > 0) samehessian = false;
			if (reset_pseudohessian) {reset_pseudohessian=false; pseudohessian = true;}

			if (accuracy < minAccuracySoFar && iterations > 0 && accuracy == std::fabs(accuracy) ) {
				minAccuracySoFar = accuracy;
			}

			if (accuracy > minAccuracySoFar*resetHessianCriterion && accuracy == std::fabs(accuracy) ) {
				if (s_info) {
					std::cout << accuracy << '\t' << minAccuracySoFar << '\t' << resetHessianCriterion << std::endl;
					std::cout << "walking backwards: newton reset" << std::endl;
				}
				resethessian(h,g,nvar);
				minAccuracySoFar *=1.5;

				if (deltamax >0.005) deltamax *=0.9;
				numIterationsSinceHessian = 0;
			}

			if (ALPHA < smallAlpha) smallAlphaCount++; else smallAlphaCount = 0;

			if (smallAlphaCount == maxNumSmallAlpha) {
				smallAlphaCount = 0;
				if (s_info) {
					std::cout << "too many small alphas: newton reset" << std::endl;
				}
				resethessian(h,g,nvar);
				if (deltamax >0.005) deltamax *=0.9;
				numIterationsSinceHessian = 0;
			}

			numIterationsSinceHessian++;
			if ((numIterationsSinceHessian >= n_iterations_for_hessian &&
						iterations > 0 && accuracy < minAccuracyForHessian && minimum < minAccuracyForHessian)) {
				if (s_info && e_info)
					std::cout << "Still no solution, computing full hessian..." << std::endl;
				else std::cout <<"*";
				pseudohessian = false; reset_pseudohessian =true;
				numIterationsSinceHessian = 0;
			}

		break;
	}
}
