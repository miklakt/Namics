#include "solve_scf.h"
#include "solve_bfgs.h"


Solve_BFGS :: Solve_BFGS(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_,vector<State*> Sta_,vector<Reaction*> Rea_,vector<Molecule*> Mol_,System* Sys_,vector<Variate*> Var_,string name_) : 
                   Solve_scf(In_,Lat_,Seg_,Sta_,Rea_,Mol_,Sys_,Var_,name_) {	

}

Solve_BFGS::~Solve_BFGS() {}


Real Solve_BFGS::operator()(Vector& x_, Vector& g_){	
	Real residual=0;
	Real*x=&x_[0];
	Real*g=&g_[0];
	residuals(x,g);
	residual=g_.norm();
       return residual;
}


//template <typename Foo>
//	Sys->CheckResults(report_errors);
//}
