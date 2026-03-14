#ifndef LGRAD2_H
#define LGRAD2_H
#include "lattice.h"
class LGrad2 : public Lattice
{
	public:	LGrad2(const Input& In_,const std::string& name_);
	~LGrad2();

	void ComputeLambdas(void);
	bool PutM();
	void TimesL(Real*);
	void DivL(Real*);
	Real Moment(Real*,Real,int);
	Real MomentPlanar(Real*,int,Real);
	Real WeightedSum(Real*);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void LReflect(Real*,Real*,Real*);
	void UReflect(Real*,Real*,Real*);
	virtual void propagateF(Real*,Real*,Real*, int, int,int);
	virtual void propagateB(Real*,Real*,Real*, int, int,int);
	bool ReadRange(int*, int*, int&, bool&, std::string, int, std::string, std::string);
	bool ReadRangeFile(std::string,int* H_p,int&, std::string, std::string);
	bool CreateMASK(Real*, int*, int*, int, bool);
	Real ComputeTheta(Real*);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_M_bounds(Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	Real ComputeGN(Real*,int, int);
	void AddPhiS(Real*,Real*,Real*,int, int);
	void AddPhiS(Real*,Real*,Real*, Real,int, int);
	void Initiate(Real*,Real*,int,int);
	void Terminate(Real*,Real*,int,int);
	void set_bounds_x(Real*,Real*,int);
	void set_bounds_y(Real*,Real*,int);
	void set_bounds_x(Real*,int);
	void set_bounds_y(Real*,int);
	bool PutMask(Real* ,std::vector<int>,std::vector<int>,std::vector<int>,int);
};
#endif
