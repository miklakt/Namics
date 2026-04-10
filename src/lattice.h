#ifndef LATTICExH
#define LATTICExH

#include "namics.h"
#include "input.h"

struct LatticeSelection {
	int gradients = 1;
	std::string geometry = "planar";
};

class Lattice {
public:
	Lattice(const Input&,const std::string&);

	virtual ~Lattice();

	std::string name;
	const Input* In;
	int MX,MY,MZ;
	std::vector<int> mx;
	std::vector<int> my;
	std::vector<int> mz;
	std::vector<int> m;
	std::vector<int> jx;
	std::vector<int> jy;
	std::vector<int> n_box;
	std::vector<Real> H;
	std::vector<std::string> BC;
	int BX1,BY1,BZ1,BXM,BYM,BZM;
	std::vector<int> B_X1;
	std::vector<int> B_Y1;
	std::vector<int> B_Z1;
	std::vector<int> B_XM;
	std::vector<int> B_YM;
	std::vector<int> B_ZM;
	int JX,JY,M;
	int JZ=1;
	bool all_lattice;
	bool ignore_sites;
	bool fcc_sites;
	bool stencil_full;
	int sub_box_on;
	int subl;
	Real volume;
	Real Accesible_volume;

	LatticeType lattice_type;
	int gradients;
	std::string geometry;
	Real offset_first_layer;
	Real bond_length;
	std::vector<Real> L;
	int Z;
	Real lambda;
	std::vector<Real> lambda0;
	std::vector<Real> fcc_lambda0;
	std::vector<Real> lambda_1;
	std::vector<Real> fcc_lambda_1;
	std::vector<Real> lambda1;
	std::vector<Real> fcc_lambda1;
	std::vector<Real> LAMBDA;
	int fjc, FJC;
	std::vector<Real> X;

	ParameterStore OUTPUT;

	void DeAllocateMemory(void);
	void AllocateMemory(void);
	void PushOutput();
	ParameterStore FormatProfile(std::span<const Real>, bool);
	std::span<Real> GetPointer(int);
	int P(int,int,int);
	int P(int,int);
	int P(int);
	bool CheckInput(int);

	void PrepareForCalculations(void);

protected:
	bool AssignChoice(const std::string&, std::string&, std::initializer_list<const char*>, const char*) const;
	bool ReadBoundaryCondition(const ParameterStore&, const char*, int, std::initializer_list<const char*>, const char*, const char* fallback = "mirror");
	bool ReadScaledDimension(const ParameterStore&, const char*, int&, int, const char*, const char*);
	void ReadOffsetFirstLayer(const ParameterStore&);
	bool RejectParameters(const ParameterStore&, std::initializer_list<std::pair<const char*, const char*>>) const;
	bool RejectAxisBoundsIn1D(const ParameterStore&) const;
	bool RejectScalarBoundsInMultiD(const ParameterStore&) const;
	bool RejectZBoundsIn2D(const ParameterStore&) const;
	bool CheckPeriodicPair(int, int, const char*) const;
	virtual bool CheckLatticeInput(const ParameterStore&) = 0;

public:
	virtual void ComputeLambdas(void) {}
	virtual Real WeightedSum(Real*)=0;
	virtual Real Moment(Real*,Real,int) =0;
	virtual void PutM(void)=0;
	virtual void propagate(Real*,Real*, int, int,int)=0;
	virtual void Side(Real *, Real *, int) =0;
	virtual void UpdateEE(Real*, Real*) =0;
	virtual void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool)=0;
	virtual void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool)=0;
	virtual void remove_bounds(Real*)=0;
	virtual void set_bounds(Real*)=0;
	virtual void remove_bounds(int*)=0;
	virtual void set_bounds(int*)=0;
	virtual void set_M_bounds(Real*)=0;
	virtual Real ComputeGN(Real*)=0;
	virtual void Initiate(Real*,Real*) =0;
};

namespace lattice_factory {
std::unique_ptr<Lattice> CreateChecked(const Input&, const std::string&, int);
}

#endif
