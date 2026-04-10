#ifndef SFNEWTONxH
#define SFNEWTONxH
#include <vector>
#include "namics.h"

class SFNewton {
public:
	SFNewton();

	virtual ~SFNewton() = default;
	bool max_g;
	int nbits;
	int lineiterations,numIterationsSinceHessian,resetiteration;
	int linesearchlimit;
	Real smallAlpha;
	int maxNumSmallAlpha;
	int smallAlphaCount;
	Real minAccuracyForHessian;
	Real minAccuracySoFar;
	Real resetHessianCriterion;
	Real trustfactor;
	bool reset_pseudohessian;
	bool pseudohessian;
	bool hessian;
	bool samehessian;
	Real accuracy;
	Real max_accuracy_for_hessian_scaling;
	int n_iterations_for_hessian;

	int trouble;
	bool e_info;
	bool hs_info;
	bool s_info;
	bool t_info;

	int i_info,iv;
	Real residual;
	std::vector<int> mask;
	int IV;
	int iterations;
	Real minimum;

	virtual void residuals(Real*,Real*) = 0; //x,g
	virtual void inneriteration(Real*,Real*, Real, Real&, Real, int) = 0; //g accuracy nvar
	void COMPUTEG(Real*,Real*,int,bool);
	void ResetX(Real*,int,bool);
	bool Message(bool,bool,int, int,Real, Real,std::string);

	Real newdirection(Real*, Real*,Real*, Real*,Real*, Real*, int, Real,bool); //there is only one of this.
	void direction(Real*, Real*, Real*, Real*, Real*, int, Real,Real,bool);
	void newhessian(Real*,Real*,Real*,Real*,Real*,int,Real,Real,bool);
	void resethessian(Real*, Real*, int);
	void startderivatives(Real*,Real*,int);
	void newtrustregion(Real*,Real,Real&,Real&,Real,Real,int); //there is only one.
	Real linesearch(Real*,Real*,Real*,Real*,Real*,int, Real,bool);  //there is only one.
	Real zero(Real*,Real*,Real*,Real*,Real*,int,Real,bool);
	Real stepchange(Real*,Real*,Real*,Real*,Real*,int,Real&,bool);
	Real linecriterion(Real*, Real*, int);
	void numhessian(Real*, Real*, Real*, int,bool);
	void findhessian(Real* ,Real*,Real*,int,bool);
	void decomposition(Real*,int,int&);
	Real norm2(Real*,int);
	void decompos(Real*, int, int&);
	int signdeterminant(Real*, int);
	void updateneg(Real* ,Real* , int, Real);
	void updatpos(Real*, Real*, Real*, int, Real);
	void gausa(Real*, Real*, Real*, int);
	void gausb(Real*, Real*, int);
	Real newfunction(Real*, int);
	Real residue(Real*, Real*, Real*, int);

	bool iterate(Real*,int,int,Real,Real,Real,bool);
	bool iterate_DIIS(Real*,int,int,int, Real, Real,int);
	void Ax(Real*, Real*, int);
	void DIIS(Real* , Real*, Real*, Real* , Real* ,Real*, int, int , int, int);
private:
	Real computeresidual(Real*, int);

};
#endif
