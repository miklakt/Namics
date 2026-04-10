#ifndef LGRAD1_H
#define LGRAD1_H
#include "lattice.h"
class LGrad1 : public Lattice
{
	public: LGrad1(const Input& In_,const std::string& name_);
	static bool Matches(const LatticeSelection& selection) {
		return selection.gradients == 1;
	}

	protected:
	bool CheckLatticeInput(const ParameterStore&) override;

	public:

	void PutM();
	Real Moment(Real*,Real,int);
	Real WeightedSum(Real*);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_bounds(Real*,Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	void set_M_bounds(Real*);

	virtual void ComputeLambdas(void);
	virtual void UpdateEE(Real*, Real*);
	virtual void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	virtual void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	virtual void Side(Real *, Real *, int);
	virtual void propagate(Real*,Real*, int, int,int);
	virtual Real ComputeGN(Real*);
	virtual void Initiate(Real*,Real*);
};
#endif
