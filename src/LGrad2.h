#ifndef LGRAD2_H
#define LGRAD2_H
#include "lattice.h"
class LGrad2 : public Lattice
{
	public:	LGrad2(const Input& In_,const std::string& name_);
	~LGrad2();
	static bool Matches(const LatticeSelection& selection) {
		return selection.gradients == 2;
	}

	protected:
	bool CheckLatticeInput(const ParameterStore&) override;

	public:
	void ComputeLambdas(void);
	bool PutM();
	Real Moment(Real*,Real,int);
	Real WeightedSum(Real*);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_M_bounds(Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	Real ComputeGN(Real*, int);
	void Initiate(Real*,Real*);
	void set_bounds_x(Real*,Real*,int);
	void set_bounds_y(Real*,Real*,int);
	void set_bounds_x(Real*,int);
	void set_bounds_y(Real*,int);
	bool PutMask(Real* ,std::vector<int>,std::vector<int>,std::vector<int>,int);
};
#endif
