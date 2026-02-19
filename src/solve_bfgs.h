#ifndef SOLVE_BFGSxH
#define SOLVE_BFGSxH
#include "solve_scf.h"


class Solve_BFGS : public Solve_scf 
{
public:

	Solve_BFGS(const Input*,Lattice*,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,System*,vector<Variate*>,string);

	~Solve_BFGS();

	

	Real operator()(Vector& , Vector& );


};
#endif
