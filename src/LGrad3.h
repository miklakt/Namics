#ifndef LGRAD3_H
#define LGRAD3_H
#include "lattice.h"
class LGrad3 : public Lattice
{
	public:	LGrad3(const Input& In_,const std::string& name_);
	~LGrad3();
	static bool Matches(const LatticeSelection& selection) {
		return selection.gradients == 3;
	}

	protected:
	bool CheckLatticeInput(const ParameterStore&) override;

	public:
	void ComputeLambdas(void);
	bool PutM();
	void TimesL(Real*);
	void DivL(Real*);
	Real Moment(Real*,Real,int);
	Real MomentPlanar(Real*,int,Real);
	Real WeightedSum(Real*);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void propagateF(Real*,Real*,Real*, int, int,int);
	void propagateB(Real*,Real*,Real*, int, int,int);
	Real ComputeTheta(Real*);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_M_bounds(Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	Real ComputeGN(Real*,int,int);
	void AddPhiS(Real*,Real*,Real*,int,int);
	void AddPhiS(Real*,Real*,Real*,Real,int,int);
	void Initiate(Real*,Real*,int,int);
	void Terminate(Real*,Real*,int,int);
	void set_bounds_x(Real*,Real*,int,int);
	void set_bounds_y(Real*,Real*,int,int);
	void set_bounds_z(Real*,Real*,int,int);
	void set_bounds_x(Real*,int,int);
	void set_bounds_y(Real*,int,int);
	void set_bounds_z(Real*,int,int);
	bool PutMask(Real* ,std::vector<int>,std::vector<int>,std::vector<int>,int);
};
#endif
