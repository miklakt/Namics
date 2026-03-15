#ifndef LGRAD1_H
#define LGRAD1_H
#include "lattice.h"
class LGrad1 : public Lattice
{
	public: LGrad1(const Input& In_,const std::string& name_);
	virtual ~LGrad1();
	static bool Matches(const LatticeSelection& selection) {
		return selection.gradients == 1 && selection.geometry != "planar";
	}


	bool PutM();
	void TimesL(Real*);
	void DivL(Real*);
	Real Moment(Real*,Real,int);
	Real MomentPlanar(Real*,int,Real);
	Real WeightedSum(Real*);
	Real ComputeTheta(Real*);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_bounds(Real*,Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	void set_M_bounds(Real*);

	virtual void ComputeLambdas(void);
	virtual void UpdateEE(Real*, Real*,Real*);
	virtual void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	virtual void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	virtual void Side(Real *, Real *, int);
	virtual void propagate(Real*,Real*, int, int,int);
	virtual void propagateF(Real*,Real*,Real*, int, int,int);
	virtual void propagateB(Real*,Real*,Real*, int, int,int);
	virtual Real ComputeGN(Real*,int,int);
	virtual void AddPhiS(Real*,Real*,Real*,int, int);
	virtual void AddPhiS(Real*,Real*,Real*,Real,int, int);
	virtual void Initiate(Real*,Real*,int, int);
	virtual void Terminate(Real*,Real*,int,int);
	bool PutMask(Real* ,std::vector<int>,std::vector<int>,std::vector<int>,int);
};
#endif
