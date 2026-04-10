#ifndef LGRAD3_H
#define LGRAD3_H
#include "lattice.h"
class LGrad3 : public Lattice
{
	public:	LGrad3(const Input& In_,const std::string& name_);
	static bool Matches(const LatticeSelection& selection) {
		return selection.gradients == 3;
	}

	protected:
	bool CheckLatticeInput(const ParameterStore&) override;

	public:
	void PutM();
	Real Moment(Real*,Real,int);
	Real WeightedSum(Real*);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void UpdateEE(Real*, Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_M_bounds(Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	Real ComputeGN(Real*);
	void Initiate(Real*,Real*);
	void set_bounds_x(Real*,Real*,int,int);
	void set_bounds_y(Real*,Real*,int,int);
	void set_bounds_z(Real*,Real*,int,int);
	void set_bounds_x(Real*,int,int);
	void set_bounds_y(Real*,int,int);
	void set_bounds_z(Real*,int,int);
};
#endif
