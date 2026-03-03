#include "solve_scf.h"
#include <iostream>

Solve_scf::Solve_scf(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_,System* Sys_,string name_) :
	name{name_}, In{In_}, Sys{Sys_}, Seg{Seg_}, Lat{Lat_}, Mol{Mol_}, Sta{Sta_}, Rea{Rea_}
{
NAMICS_DBG("Constructor in Solve_scf " << endl);
	lat=Lat;
	KEYS.push_back("gradient_type");
	KEYS.push_back("method");
	KEYS.push_back("x_info");
	KEYS.push_back("e_info"); KEYS.push_back("s_info");KEYS.push_back("i_info");KEYS.push_back("t_info");KEYS.push_back("hs_info");
	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance");
	KEYS.push_back("stop_criterion");
	KEYS.push_back("deltamin");KEYS.push_back("deltamax");
	KEYS.push_back("linesearchlimit");
	//KEYS.push_back("samehessian");
	KEYS.push_back("max_accuracy_for_hessian_scaling");
	KEYS.push_back("n_iterations_for_hessian");
	KEYS.push_back("small_alpha");
	//KEYS.push_back("target_function");
	KEYS.push_back("max_n_small_alpha");
	KEYS.push_back("min_accuracy_for_hessian");
	KEYS.push_back("max_fr_reverse_direction");
	KEYS.push_back("print_hessian_at_it");
	KEYS.push_back("m");
	KEYS.push_back("n_restart_DIIS");
	max_g = false; // compute g based on max error
	rescue_status = NONE;
	all=false;
	restart_DIIS =0;
}

Solve_scf::~Solve_scf() {
	DeAllocateMemory();
}

void Solve_scf :: DeAllocateMemory(){
NAMICS_DBG("DeAllocateMemory in Solve " << endl);

		int niv = In->ReactionList.size();
		if (niv>0) {
			free(yy);
			free(SIGN);
		}
		//delete [] xx;
		//free(xx);
all=false;
NAMICS_DBG("exit for 'destructor' in Solve " << endl);

}

void Solve_scf::AllocateMemory() {
NAMICS_DBG("AllocateMemeory in Solve " << endl);
	if (all) DeAllocateMemory();
	int M=lat->M;
	iv = (Sys->ItMonList.size() + Sys->ItStateList.size())* M;
	if (Sys->charged) iv += M;
	if (SCF_method=="Picard") iv += M;
	if (Sys->constraintfields) iv +=M;
	int length = In->MonList.size();
	for (int i = 0; i < length; i++) iv+=Seg[i]->constraint_z.size();
	x_storage.assign(iv, 0);
	xx=x_storage.data();
	all=true;
	int niv = In->ReactionList.size();
	if (niv>0) {
		yy=(Real*) malloc(niv*sizeof(Real)); std::fill_n(yy, niv, 0);
		SIGN=(int*) malloc((niv)*sizeof(int)); for (int i=0; i<niv; i++) SIGN[i]=1.0;
	}

	Sys->AllocateMemory();
}

bool Solve_scf::PrepareForCalculations() {
NAMICS_DBG("PrepareForCalculations in Solve " << endl);
	bool success=true;
	return success;
}

bool Solve_scf::CheckInput(int start_) { start=start_;
NAMICS_DBG("CheckInput in Solve " << endl);
	pseudohessian =false;
	deltamin =0.1;
	s_info=false;
	e_info=false;
	t_info=false;
	hs_info=false;
	i_info=1;
	hessian =false;
	bool success=true;
	control=proceed;
	string value;
	solver=PSEUDOHESSIAN;
	SCF_method="pseudohessian";
	gradient=classical;
	residual=1;
	m=10;
	success=In->CheckParameters("newton",name,start, KEYS, PARAMETERS);
	if (success) {
		iterationlimit=ParseInt(GetValue("iterationlimit"),1000);
		if (iterationlimit < 0 || iterationlimit>1e6) {iterationlimit = 1000;}

		e_info=ParseBool(GetValue("e_info"),true); value_e_info=e_info;
		hs_info=ParseBool(GetValue("hs_info"),true);
		s_info=ParseBool(GetValue("s_info"),false); value_s_info =s_info;
		t_info=ParseBool(GetValue("t_info"),false);
		i_info=ParseInt(GetValue("i_info"),1);
		if (i_info == 0) {
		// We cannot divide by zero (see modulus statements in sfnewton), but this will probably be what the user means.
		cerr << "WARNING: i_info cannot be zero ! Defaulting to iterationlimit + 1."<< endl;
		i_info = iterationlimit+1;
		}
		value_i_info=i_info;

		if (GetValue("target_function").size() > 0) {
			string target;
      			target = ParseString(GetValue("target_function"), target);
			using namespace std::placeholders;
			if ( target.find("log") != string::npos  ) target_function = bind(&Solve_scf::gradient_log, this, _1, _2, _3, _4, _5);
			else if ( target.find("quotient") != string::npos  ) target_function = bind(&Solve_scf::gradient_quotient, this, _1, _2, _3, _4, _5);
			else if ( target.find("minus") != string::npos  ) target_function = bind(&Solve_scf::gradient_minus, this, _1, _2, _3, _4, _5);
			else {
				cerr << "Target function not found, please choose from log, quotient or minus. Defaulting to minus." << endl;
			}
		} else {
			using namespace std::placeholders;
			target_function = bind(&Solve_scf::gradient_minus, this, _1, _2, _3, _4, _5);
		}

		deltamax=ParseReal(GetValue("deltamax"),0.1);
		if (deltamax < 0 || deltamax>100) {deltamax = 0.1;  cout << "Value of deltamax out of range 0..100, and value set to default value 0.1" <<endl; }
		deltamin=0;
		deltamin=ParseReal(GetValue("deltamin"),deltamin);
		if (deltamin < 0 || deltamin>100) {deltamin = deltamax/100000;  cout << "Value of deltamin out of range 0..100, and value set to default value deltamax/100000" <<endl; }
		tolerance=ParseReal(GetValue("tolerance"),1e-7);
		if (tolerance < 1e-16 ||tolerance>10) {tolerance = 1e-5;  cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<endl; }

			if (GetValue("method").size()==0) {SCF_method="pseudohessian";} else {
				vector<string>method_options;
				method_options.push_back("DIIS");
				//method_options.push_back("Picard"); //can be included again when adjusted for charges and guess
				method_options.push_back("pseudohessian");
				method_options.push_back("hessian");
				//method_options.push_back("conjugate_gradient");
				method_options.push_back("LBFGS");
				method_options.push_back("BRR");
				if (!ParseString(GetValue("method"),SCF_method,method_options,"In 'solve_scf' the entry for 'method' not recognized: choose from:")) success=false;
			}
		if (SCF_method=="hessian" || SCF_method=="pseudohessian") {
			if (SCF_method=="hessian") {pseudohessian=false; hessian=true; solver=HESSIAN;} else { pseudohessian=true; hessian=false; solver=PSEUDOHESSIAN;}
			samehessian=false; //ParseBool(GetValue("samehessian"),false);
			max_accuracy_for_hessian_scaling=ParseReal(GetValue("max_accuracy_for_hessian_scaling"),0.1);
			if (max_accuracy_for_hessian_scaling<1e-7 || max_accuracy_for_hessian_scaling>1) {
				cout <<"max_accuracy_for_hessian_scaling is out of range: 1e-7...1; default value 0.1 is used instead" << endl;
				max_accuracy_for_hessian_scaling=0.1;
			}
			minAccuracyForHessian=ParseReal(GetValue("min_accuracy_for_hessian"),0.5);
			if (minAccuracyForHessian<0 ||minAccuracyForHessian>1) {
				cout <<"min_accuracy_for_hessian is out of range: 0...0.1; default value 0 is used instead (no hessian computation)" << endl;
				minAccuracyForHessian=0;
			}
			maxFrReverseDirection =ParseReal(GetValue("max_fr_reverse_direction"),0.4);
			if (maxFrReverseDirection <0.1 ||maxFrReverseDirection >0.5) {
				cout <<"max_fr_reverse_direction is out of range: 0.1...0.5; default value 0.4 is used instead" << endl;
				maxFrReverseDirection =0.4;
			}

			n_iterations_for_hessian=ParseInt(GetValue("n_iterations_for_hessian"),iterationlimit+100);
			if (n_iterations_for_hessian<1 ) {
				cout <<" n_iterations_for_hessian setting must be larger than unity; hessian evaluations will not be done " << endl;
				n_iterations_for_hessian=iterationlimit+100;
			}
			maxNumSmallAlpha=ParseInt(GetValue("max_n_small_alpha"),50);
			if (maxNumSmallAlpha<10 ||maxNumSmallAlpha>1000) {
				cout <<" max_n_small_alpha is out of range: 10, ..., 100;  max_n_small_alpha is set to default: 50 " << endl;
				maxNumSmallAlpha=50;
			}

			deltamin=ParseReal(GetValue("delta_min"),0);
			if (deltamin <0 || deltamin>deltamax) {
				cout <<"delta_min is out of range; 0, ..., " << deltamax << "; delta_min value set to 0 " << endl;
				deltamin=0;
			}
			smallAlpha=ParseReal(GetValue("small_alpha"),0.00001);
			if (smallAlpha <0 || smallAlpha>1) {
				cout <<"small_alpha is out of range; 0, ..., 1; small_alpha value set to default: 1e-5 " << endl;
				smallAlpha=0.00001;
			}
		}
		if (SCF_method=="DIIS") {
			solver=diis;
			m=ParseInt(GetValue("m"),10);
			if (m < 0 ||m>100) {m=10;  cout << "Value of 'm' out of range 0..100, value set to default value 10" <<endl; }
			restart_DIIS=iterationlimit;
			restart_DIIS=ParseInt(GetValue("n_restart_DIIS"),iterationlimit);
			if (restart_DIIS < 0 || restart_DIIS > iterationlimit*10) {
				restart_DIIS=iterationlimit; cout <<"Value of 'n_restart_DIIS' out of range 0 .. iterationlimit; value set to iterationlimit" << endl;
			}
			restart_DIIS -=restart_DIIS%m; cout <<"Restart DIIS set to " << restart_DIIS << endl;
		}
		if (SCF_method=="Picard") {
			solver= PICARD;
			gradient=Picard;
		}
		if (SCF_method=="conjugate_gradient") {
			solver= conjugate_gradient;
			linesearchlimit=ParseInt(GetValue("linesearchlimit"),linesearchlimit);
		}

			if (SCF_method=="LBFGS") {
				solver=LBFGS;
				m=ParseInt(GetValue("m"),6);
				if (m < 0 ||m>1000) {m=6;  cout << "Value of 'm' out of range 0..1000, value set to default value 6" <<endl; }
			}

		if (SCF_method=="BRR") {
			solver=BRR;
			m=ParseInt(GetValue("m"),10);
			if (m < 0 ||m>1000) {m=10;  cout << "In method 'BRR', value of 'm' out of range 0..1000, value set to default value 10" <<endl; }
		}

			if (GetValue("gradient_type").size()==0) {gradient=classical;} else {
				vector<string>gradient_options;
				gradient_options.push_back("classical");
				//gradient_options.push_back("Picard");
				if (!ParseString(GetValue("gradient_type"),gradients,gradient_options,"In 'solve_scf' the entry for 'gradient_type' not recognized: choose from:")) success=false;
				if (gradients=="classical") gradient=classical;
				if (gradients=="Picard")  gradient=Picard;
			}

		StoreFileGuess=ParseString(GetValue("store_guess"),"");
		ReadFileGuess=ParseString(GetValue("read_guess"),"");
		if (GetValue("stop_criterion").size() > 0) {
			vector<string>options;
			options.push_back("norm_of_g");
			options.push_back("max_of_element_of_|g|");
			if (!ParseString(GetValue("stop_criterion"),stop_criterion,options,"In newton the stop_criterion setting was not recognised")) {success=false; };
			if(GetValue("stop_criterion") == options[1]) {
				max_g = true;
			}
		}
	}
	return success;
}


void Solve_scf::PutParameter(string new_param) {
NAMICS_DBG("PutParameter in Solve " << endl);
	KEYS.push_back(new_param);
}

string Solve_scf::GetValue(string parameter){
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

void Solve_scf::push(string s, Real X) {
NAMICS_DBG("push (Real) in  Solve " << endl);
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Solve_scf::push(string s, int X) {
NAMICS_DBG("push (int) in  Solve " << endl);
	ints.push_back(s);
	ints_value.push_back(X);
}
void Solve_scf::push(string s, bool X) {
NAMICS_DBG("push (bool) in  Solve " << endl);
	bools.push_back(s);
	bools_value.push_back(X);
}
void Solve_scf::push(string s, string X) {
NAMICS_DBG("push (string) in  Solve " << endl);
	strings.push_back(s);
	strings_value.push_back(X);
}
void Solve_scf::PushOutput() {
NAMICS_DBG("PushOutput in  Solve " << endl);
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("method",SCF_method);
	push("m",m);
	push("delta_max",deltamax);
	push("residual",residual);
	push("tolerance",tolerance);
	push("iterations",iterations);
	push("iterationlimit",iterationlimit);
	push("stop_criterion",stop_criterion);
	if (pseudohessian || hessian) {
		//push("same_hessian",samehessian);
		push("linesearchlimit",linesearchlimit);
		push("max_accuracy_for_hessian_scaling",max_accuracy_for_hessian_scaling);
		push("n_iteratons_for_hessian",n_iterations_for_hessian);
		push("small_alpha",smallAlpha);
		push("max_n_small_alpha",maxNumSmallAlpha);
		push("min_accuracy_for_hessian",minAccuracyForHessian);
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

int Solve_scf::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
NAMICS_DBG("GetValue (long) in  Solve " << endl);
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

void Solve_scf::Copy(Real* x, Real* X, int MX, int MY, int MZ, int fjc_old) {
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
					cout <<" Copy from more than one gradient to one gradient: (i) =(1,i) or (1,1,i) is used "<< endl;
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
				cout <<" Copy from one-gradient to two gradients: one-gradient (x,i)=(i) is used for all x " << endl;
				for (i=0; i<mx+2*fjc; i++)
				for (j=0; j<my+2*fjc; j++) if (j<MX+2*fjc_old) x[i*jx+j]=X[j];
			} else {
				if (MZ>0) {
					cout <<" Copy from three gradients to two gradients: (i,j)=(1,i,j) is used " <<endl;
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
					cout <<"Copy from one gradient to three gradients: (x,y,i) = (i) is used for all x,y " << endl;
					for (i=0; i<mx+2*fjc; i++)
					for (j=0; j<my+2*fjc; j++)
					for (k=0; k<mz+2*fjc; k++) if (k<MX+2*fjc_old) x[i*jx+j*jy+k] = X[k];
				} else {
					if (MZ==0) {
						cout <<"Copy form two gradients to three: (x,i,j) = (i,j) for all x " << endl;
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

bool Solve_scf::Guess(Real *X, string METHOD, vector<string> MONLIST, vector<string> STATELIST, bool CHARGED, int MX, int MY, int MZ,int fjc_old){
	(void)METHOD;
	NAMICS_DBG( "Guess in Solve" << endl);
	int M=lat->M;
	bool success=true;
	if (start ==1 && Sys->GuessType != "")  {
		cout <<"guessing " << endl;
		lat->GenerateGuess(xx,Sys->CalculationType,Sys->GuessType,Seg[Sys->MonA]->guess_u,Seg[Sys->MonB]->guess_u);
	} else {
		int m;
		if (MZ>0) {m=(MX+2)*(MY+2)*(MZ+2); } else { if (MY>0) { m=(MX+2*fjc_old)*(MY+2*fjc_old); } else {  m=(MX+2*fjc_old);}}

		int length_old_mon=MONLIST.size();
		int length_old_state=STATELIST.size();
		int length_new_mon=Sys->ItMonList.size();
		int length_new_state=Sys->ItStateList.size();
		for (int i = 0; i<length_old_mon; i++) {
			for (int j=0; j<length_new_mon; j++) {
				if (MONLIST[i]==Seg[Sys->ItMonList[j]]->name) {
					Copy(xx+M*j,X+i*m,MX,MY,MZ,fjc_old);
				}
			}
		}
		for (int i = 0; i<length_old_state; i++) {
			for (int j=0; j<length_new_state; j++) {
				if (STATELIST[i]==Sta[Sys->ItStateList[j]]->name) {
					Copy(xx+M*(j+length_new_mon),X+(i+length_old_mon)*m,MX,MY,MZ,fjc_old);
				}
			}
		}

		if (CHARGED && Sys->charged) {
			Copy(xx+(length_new_mon+length_new_state)*M,X+(length_old_mon+length_old_state)*m,MX,MY,MZ,fjc_old);
		}
	}
	return success;
}

class SCF_LBFGS
{
private:
    const Input* In;
    Lattice* Lat;
    vector<Segment*> Seg;
    vector<State*> Sta;
    vector<Reaction*> Rea;
    vector<Molecule*> Mol;
    System* Sys;
    int iterations =0;
    Real residual=1;
public:
    SCF_LBFGS(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_,vector<State*> Sta_,vector<Reaction*> Rea_,vector<Molecule*> Mol_,System* Sys_) :
      In(In_),Lat(Lat_),Seg(Seg_),Sta(Sta_),Rea(Rea_),Mol(Mol_),Sys(Sys_) {}

    Real operator()(Vector& x_, Vector& g_)
    {
	Real* x=&x_[0];
	Real* g=&g_[0];
	int iv=x_.size();
	Sys->Classical_residual(x,g,residual,iterations, iv);
	iterations++;
	residual=g_.norm();
	return residual;
	    }
};

bool Solve_scf::Solve(bool report_errors_) { //going SCF here
NAMICS_DBG("Solve in  Solve_scf " << endl);
	bool success=true;
	bool report_errors=report_errors_;
	int niv = In->ReactionList.size();
	if (niv>0) {
		int i_solver=0;
		if (solver==HESSIAN) i_solver=1;
			if (solver==PSEUDOHESSIAN) i_solver=2;
			if (solver==diis) i_solver=3;
			if (solver==BRR) i_solver=4;
			if (solver==LBFGS) i_solver=5;
		bool ee_info, ss_info;
		if (e_info) ee_info=true; else ee_info=false; e_info=false;
		if (s_info) ss_info=true; else ss_info=false; s_info=false;
		gradient = WEAK;
		control= super;
		pseudohessian=false; hessian =true;

		success=iterate(yy,niv,100,1e-8,1,0.0000001,true);
		cout << iterations << " iterations to find alphabulk values. " <<endl;
		if (!success) cout <<"iteration for alphabulk values for internal states failed. Check eqns. " << endl;
		e_info=ee_info;
		s_info=ss_info;
		if (i_solver==1) solver=HESSIAN;
			if (i_solver==2) {solver=PSEUDOHESSIAN; pseudohessian=true;}
			if (i_solver==3) solver=diis;
			if (i_solver==4) solver=BRR;
			if (i_solver==5) solver=LBFGS;
			gradient = classical;
			control = proceed;
		}

		switch(solver) {
		case HESSIAN:
			success=iterate(xx,iv,iterationlimit,tolerance,deltamax,deltamin,true);
		break;
		case PSEUDOHESSIAN:
			success=iterate(xx,iv,iterationlimit,tolerance,deltamax,deltamin,true);
		break;
		case PICARD:
			success=iterate_Picard(xx,iv,iterationlimit,tolerance,deltamax);
		break;
		case diis:
			success=iterate_DIIS(xx,iv,m,iterationlimit,tolerance,deltamax,restart_DIIS);
		break;
		case BRR:
			success=iterate_BRR(xx,iv,m,iterationlimit,tolerance,deltamax);
		break;
		case conjugate_gradient:
			success =iterate_conjugate_gradient(xx,iv,iterationlimit,tolerance,deltamax);
		break;
		case LBFGS:
			success=true;
			{
			SCF_LBFGS fun(In,Lat,Seg,Sta,Rea,Mol,Sys);
			LBFGSParam<Real> param;
			param.epsilon=tolerance;
			param.m=m;
			param.max_iterations =iterationlimit;
			if (deltamax > 0) param.max_step = deltamax;
			if (deltamin > 0 && deltamax > deltamin) param.min_step = deltamin;
			LBFGSSolver<Real> mysolver(param);
			Real fx=0;
			cout <<endl <<"LBFGS has been notified" << endl;
			Vector x_vec = Eigen::Map<Vector>(xx, iv);
			iterations =mysolver.minimize(fun, x_vec, fx);
			std::copy_n(x_vec.data(), iv, xx);
			Real res = mysolver.final_grad_norm();
			cout <<endl <<"Problem solved: " << iterations << " iterations,  |g|: " << res <<  endl;
			}
		break;
		default:
			cout <<"Solve is lost" << endl; success=false;
		break;
	}
	success=Sys->CheckResults(report_errors);
	return success;
}


void Solve_scf::residuals(Real* x, Real* g){
 NAMICS_DBG("residuals in Solve_scf " << endl);
	int M=lat->M;
	Real chi;
	int sysmon_length = Sys->SysMonList.size();
	int mon_length = In->MonList.size(); //also frozen segments

	switch(gradient) {
		case WEAK:
			NAMICS_DBG("Residuals for weak iteration " << endl);
			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				if (Rea[i]->Sto.size()==3) Rea[i]->GuessAlpha();
			}
			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				if (Rea[i]->Sto.size()!=3) Rea[i]->GuessAlpha();
			}
			for (size_t i = 0; i<In->ReactionList.size(); i++) {
				Rea[i]->PutAlpha(exp(x[i]));
			}


			std::fill_n(g, In->ReactionList.size(), 0);

			for (size_t i = 0; i<In->ReactionList.size(); i++) {

				g[i]=SIGN[i]*Rea[i]->Residual_value();
				//g[i]=Rea[i]->Residual_value();

			}

		break;
		case Picard:
		{
			NAMICS_DBG("Residuals in Picard mode in Solve_scf " << endl);
			int jump=sysmon_length;
			if (Sys->charged) jump++;
			std::copy_n(xx+jump*M, M, alpha);
			Sys->ComputePhis(x,iterations==0,residual);
			if (Sys->charged) {
				Sys->DoElectrostatics(g+sysmon_length*M,xx+sysmon_length*M);
				lat->UpdateEE(Sys->EE,Sys->psi,Sys->E);
				lat->set_bounds(Sys->psi);
				lat->UpdatePsi(g+sysmon_length*M,Sys->psi,Sys->q,Sys->eps,Sys->psiMask,Sys->grad_epsilon,Sys->fixedPsi0);
				lat->remove_bounds(g+sysmon_length*M);
			}
			Real one=1.0;
			for (int __i = 0; __i < (M); ++__i) (g+jump*M)[__i] = (Sys->phitot)[__i] + (-1.0*one);
			for (int i=0; i<sysmon_length; i++) {
				std::copy_n(xx+i*M, M, g+i*M);
				for (int k=0; k<mon_length; k++) {
                       		chi= -1.0*Sys->CHI[Sys->SysMonList[i]*mon_length+k];  //The minus sign here is to change the sign of x! just a trick due to properties of PutAlpha where a minus sing is implemented....
					if (chi!=0) for (int __i = 0; __i < (M); ++__i) if ((Sys->phitot)[__i] > 0) (g+i*M)[__i] = (g+i*M)[__i] - (chi) * (((Seg[k]->phi_side)[__i] / (Sys->phitot)[__i]) - (Seg[k]->phibulk));
				}
				if (Sys->charged){
					for (int __i = 0; __i < (M); ++__i) (g+i*M)[__i] += (Seg[Sys->SysMonList[i]]->epsilon) * (Sys->EE)[__i];
					if (Seg[Sys->SysMonList[i]]->valence !=0)
					for (int __i = 0; __i < (M); ++__i) (g+i*M)[__i] += (-1.0*Seg[Sys->SysMonList[i]]->valence) * (Sys->psi)[__i];
				}
				lat->remove_bounds(g+i*M);
				for (int __i = 0; __i < (M); ++__i) (g+i*M)[__i] = (g+i*M)[__i] * (Sys->KSAM)[__i];
			}
		break;
		}
		default:
			NAMICS_DBG("Residuals in scf mode in Solve_scf " << endl);
			Sys->Classical_residual(x,g,residual,iterations, iv);
		break;
	}
}

void Solve_scf::gradient_log(Real* g, int k, int M, int i, int j) {
	//Target function: ln(g/phi) < tolerance
	for (int z = 0 ; z < M ; ++z) {
		Real frac = (g+k*M)[z]/(Mol[i]->phi+j*M)[z];
		(g+k*M)[z] = log( frac );
	}
}

void Solve_scf::gradient_quotient(Real* g, int k, int M, int i, int j) {
	//Target function: g/phi-1 < tolerance
	for (int z = 0 ; z < M ; ++z) {
		Real frac = (g+k*M)[z]/(Mol[i]->phi+j*M)[z];
		(g+k*M)[z] = frac - 1;
	}
}

void Solve_scf::gradient_minus(Real* g, int k, int M, int i, int j) {
	//Target function: g - phi < tolerance
		for (int __i = 0; __i < (M); ++__i) (g+k*M)[__i] -= (Mol[i]->phi+j*M)[__i];
}

void Solve_scf::inneriteration(Real* x, Real* g, Real* h, Real accuracy, Real& deltamax, Real ALPHA, int nvar) {
NAMICS_DBG("inneriteration in Solve_scf " << endl);
	residual=accuracy; //hoping this is not creating problems with the use of residual...
	switch(control) {
		case super:
			for (int i=0; i<nvar; i++) { //this is to control the sing in the WEAK iteration.
				if (h[i+i*nvar]<0) {SIGN[i] = -1;  }
			}

		break;
		default:
			if (iterations > 0) samehessian = false;
			if (reset_pseudohessian) {reset_pseudohessian=false; pseudohessian = true;}

			if (accuracy < minAccuracySoFar && iterations > 0 && accuracy == fabs(accuracy) ) {
				minAccuracySoFar = accuracy;
			}

			if (accuracy > minAccuracySoFar*resetHessianCriterion && accuracy == fabs(accuracy) ) {
				if (s_info) {
					cout << accuracy << '\t' << minAccuracySoFar << '\t' << resetHessianCriterion << endl;
					cout << "walking backwards: newton reset" << endl;
				}
				resethessian(h,g,x,nvar);
				minAccuracySoFar *=1.5;

				if (deltamax >0.005) deltamax *=0.9;
				numIterationsSinceHessian = 0;
			}

			if (ALPHA < smallAlpha) smallAlphaCount++; else smallAlphaCount = 0;

			if (smallAlphaCount == maxNumSmallAlpha) {
				smallAlphaCount = 0;
				if (s_info) {
					cout << "too many small alphas: newton reset" << endl;
				}
				resethessian(h,g,x,nvar);
				if (deltamax >0.005) deltamax *=0.9;
				numIterationsSinceHessian = 0;
			}

			if (!newtondirection && pseudohessian) {
				reverseDirection[iterations%reverseDirectionRange] = 1;
			} else {
				reverseDirection[iterations%reverseDirectionRange] = 0;
			}

			numReverseDirection = 0;
			for (int i=0; i<reverseDirectionRange; i++) {
				if (reverseDirection[i] == 1)
					numReverseDirection++;
			}

			numIterationsSinceHessian++;
			Real frReverseDirection = Real(numReverseDirection)/reverseDirectionRange;
			if ((frReverseDirection > maxFrReverseDirection && pseudohessian && accuracy < minAccuracyForHessian)) {
				if (s_info && e_info) cout <<"Bad convergence (reverse direction), computing full hessian..." << endl; else cout <<"!";
				pseudohessian = false; reset_pseudohessian =true;
				numIterationsSinceHessian = 0;
			} else if ((numIterationsSinceHessian >= n_iterations_for_hessian &&
						iterations > 0 && accuracy < minAccuracyForHessian && minimum < minAccuracyForHessian)) {
				if (s_info && e_info)
					cout << "Still no solution, computing full hessian..." << endl;
				else cout <<"*";
				pseudohessian = false; reset_pseudohessian =true;
				numIterationsSinceHessian = 0;
			}

		break;
	}
}

bool Solve_scf::attempt_DIIS_rescue() {
	cout << "Attempting rescue!" << endl;
	switch (rescue_status) {
		case NONE:
			cout << "Zeroing iteration variables." << endl;
			std::fill_n(xx, iv, 0);
			rescue_status = ZERO;
			break;
		case ZERO:
			cout << "Adjusting memory depth." << endl;
			m *= 0.5;
			cout << "Zeroing iteration variables." << endl;
			std::fill_n(xx, iv, 0);
			rescue_status = M;
			break;
		case M:
			cout << "Decreasing delta_max." << endl;
			deltamax *= 0.1;
			cout << "Zeroing iteration variables." << endl;
			std::fill_n(xx, iv, 0);
			rescue_status = DELTA_MAX;
			break;
		case DELTA_MAX:
			cerr << "Exhausted all rescue options. Crash is imminent, exiting." << endl;
			exit(0);
			break;
	}
	return true;
}
