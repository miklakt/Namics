#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include "sfnewton.h"
#include "tools.h"


SFNewton::SFNewton () : residual{0} {

/*      class for
        unconstrained minimization,
        systems of nonlinear algebraic equations,
        curve fitting.

References:

 (1)    "The State of the Art in Numerical Analysis"
        Edited by D.Jacobs. (1977) Academic Press.
 (2)    "Numerical Methods for Unconstrained Optimization"
        Edited by W.Murray. (1972) Academic Press.
 (3)    "Numerical Methods for Constrained Optimization"
        Edited by P.E.Gill and W.Murray. (1974) Academic Press.
 (4)    "Numerical Methods for Non-linear Algebraic Equations"
        Edited by P.Rabinowitz. (1970) Gordon and Breach.
 (5)    "Minimization Subject to Bounds on the Variables"
        P.E.Gill and W.Murray. NPL Report NAC 72 (1976).
 (6)    "Safeguarded Steplength Algorithms for Optimization Using
        Descent Methods" P.E.Gill and W.Murray. NPL Report NAC 37
        (1974).
 (7)    "The Implementation of Two Modified Newton Methods for
        Unconstrained Optimization"
        P.E.Gill, W.Murray, and S.M.Picken. NPL Rpt. NAC 24  (1972)
 (8)    "The Implementation of Two Quasi-Newton Methods for
        Unconstrained Optimization"
        P.E.Gill, W.Murray, and R.A.Pitfield. NPL Rpt NAC 11 (1972)
 (9)    "Conjugate Gradient Methods with Inexact Searches"
        D.F.Shanno. Math. of Operations Res. 3 (1978) 244-256

Author:
Jan Scheutjens (1947-1992), Wageningen Agricultural University, NL.

C Copyright (1980) (1981-1989) Wageningen Agricultural University, NL.

C++ translation:
Peter Barneveld, Wageningen Agricultural University, NL.
Addaptation for using std::vector class in namics (filtering with mask):
Frans Leermakers, Wageningen Agricultural University, NL.

C Copyright (2018) Wageningen University, NL.

 *NO PART OF THIS WORK MAY BE REPRODUCED, EITHER ELECTRONICALLY OF OTHERWISE*

*/
       nbits = std::numeric_limits<Real>::digits;	//nbits=52;
       linesearchlimit = iterations=lineiterations=numIterationsSinceHessian = trouble=resetiteration=0;
       trouble = resetiteration = 0;
	numReverseDirection =0;
	trustregion =0.0;
	pseudohessian = samehessian = false;
	e_info = s_info = false;
	newtondirection  = false ;
	ignore_newton_direction = true;
	i_info=1;
	max_accuracy_for_hessian_scaling = 0.1;
	linesearchlimit = 20;
	linetolerance = 9e-1;
	epsilon = 0.1/std::pow(2.0,nbits/2);
	minAccuracySoFar = 1e30;
	reverseDirectionRange = 50;
	resetHessianCriterion = 1e5;
	reset_pseudohessian = false;
	accuracy=1e30;
	numIterationsSinceHessian=0;
	maxFrReverseDirection=0.4;
	numIterationsForHessian=100;
	minAccuracyForHessian=0.1;
	reverseDirection.assign(reverseDirectionRange, 0);

}

SFNewton::~SFNewton() {
}


void SFNewton::multiply(Real *v,Real alpha, Real *h, Real *w, int nvar) { //done
NAMICS_DBG("multiply in Newton" << std::endl);
	int i=0,i1=0,j=0;
	Real sum=0;
	std::vector<Real> x(nvar, 0);
	for (i=0; i<nvar; i++) {
		sum = 0;
		i1 = i-1;
		for (j=i+1; j<nvar; j++) {
			sum += w[j] * h[i+nvar*j];
		}
		x[i] = (sum+w[i])*h[i+nvar*i];
		sum = 0;
		for (j=0; j<=i1; j++) {
			sum += x[j] * h[i+nvar*j];
		}
		v[i] = alpha*(sum+x[i]);
	}
}


Real SFNewton::norm2(Real*x, int nvar) { //done
NAMICS_DBG("norm2 in Newton" << std::endl);
	const Real sum = std::inner_product(x, x + nvar, x, Real(0));
	return std::sqrt(sum);
}

int SFNewton::signdeterminant(Real*h,int nvar) { //dome
NAMICS_DBG("signdeterminant in Newton" << std::endl);
	int sign=1;
	for (int i=0; i<nvar; i++) {
		if ( h[i+i*nvar]<0 ) {
			sign = -sign;
		}
	}
	return sign;
}

void SFNewton::updateneg(Real *l,Real *w, int nvar, Real alpha) { //done
NAMICS_DBG("updateneg in Newton" << std::endl);
	int i=0,i1=0,j=0;
	Real dmin=0,sum=0,b=0,d=0,p=0,lji=0,t=0;
	dmin = 1.0/std::pow(2.0,54);
	alpha = std::sqrt(-alpha);
	for (i=0; i<nvar; i++) {
		i1 = i-1;
		sum = 0;
		for (j=0;j<=i1; j++) {
			sum += l[i+nvar*j]*w[j];
		}
		w[i] = alpha*w[i]-sum;
		t += (w[i]/l[i+nvar*i])*w[i];
	}
	t = 1-t;
	if ( t<dmin ) t = dmin;
	for (i=nvar-1; i>=0; i--) {
		p = w[i];
		d = l[i+nvar*i];
		b = d*t;
		t += (p/d)*p;
		l[i+nvar*i] = b/t;
		b = -p/b;
		for (j=i+1; j<nvar; j++) {
			lji = l[j+nvar*i];
			l[j+nvar*i] = lji+b*w[j];
			w[j] += p*lji;
		}
	}
}

void SFNewton::decompos(Real *h, int nvar, int &ntr) { //done
NAMICS_DBG("decompos in Newton" << std::endl);
	int i,j,k;//itr,ntr;
	Real sum,lsum,usum,phi,phitr,c,l;
	Real *ha,*hai,*haj;
	ha = &h[-1];
	phitr = std::numeric_limits<Real>::max();
	ntr = 0;
	i = 0;
	while (i++<nvar) {
		hai = &ha[(i-1)*nvar];
		sum = 0;
		j = 0;
		while (j++<i-1) {
			haj = &ha[(j-1)*nvar];
			c = haj[i];
			l = c/haj[j];
			haj[i] = l;
			c = hai[j];
			hai[j] = c/haj[j];
			sum += l*c;
		}
		phi = hai[i] - sum;
		hai[i] = phi;
		if (phi<0) ntr++;
		if (phi<phitr) phitr = phi;
		j = i;
		while (j++<nvar) {
			haj = &ha[(j-1)*nvar];
			lsum = 0;
			usum = 0;
			k = 0;
			while (k++<i-1) {
				lsum += ha[(k-1)*nvar+j]*hai[k];
				usum += ha[(k-1)*nvar+i]*haj[k];
			}
			hai[j] -= lsum;
			haj[i] -= usum;
		}
	}
}

void SFNewton::updatpos(Real *l, Real *w, Real *v, int nvar, Real alpha) { //done
NAMICS_DBG("updatepos in Newton" << std::endl);
	int i,j;
	Real b,c,d;
	Real vai,waj,vaj;
	Real *lai,*laj;
	Real * wa = &w[-1];
	Real * va = &v[-1];
	i = 0;
	while (i++<nvar) {
		vai = va[i];
		lai = &l[-1 + (i-1)*nvar];
		d = lai[i];
		b = d+(alpha*wa[i])*vai;
		lai[i] = b;
		d /= b;
		c = vai*alpha/b;
		b = wa[i]*alpha/b;
		alpha *= d;
		j = i;
		while (j++<nvar) {
			waj = wa[j];
			wa[j] -= wa[i]*lai[j];
			lai[j] *= d;
			lai[j] += c*waj;
			laj = &l[-1 + (j-1)*nvar];
			vaj = va[j];
			va[j] -= vai*laj[i];
			laj[i] *= d;
			laj[i] += b*vaj;
		}
	}
}

void SFNewton::gausa(Real *l, Real *dup, Real *g, int nvar) {//done
NAMICS_DBG("gausa in Newton" << std::endl);
	int i,j;
	Real*dupa,sum;
	Real *ga;
	Real *lai;

	dupa = &dup[-1];
	ga = &g[-1];

	i = 0;
	while (i++<nvar) {
		sum = 0;
		lai = &l[i-1];
		j = 0;
		while (j++<i-1) {
			sum += lai[(j-1)*nvar]*dupa[j];
		}
		dupa[i] = - ga[i] - sum;
	}
}

void SFNewton::gausb(Real *du, Real *p, int nvar) { //done
NAMICS_DBG("gausb in Newton " << std::endl);
	int i,j;
	Real *pa,sum;
	Real *duai;
	pa = &p[-1];
	i = nvar+1;
	while (i-- > 1) {
		sum = 0;
		duai = &du[i-1];
		j = i;
		while (j++<nvar) {
			sum += duai[(j-1)*nvar]*pa[j];
		}
		pa[i] = pa[i]/duai[(i-1)*nvar] - sum;
	}
}

Real SFNewton::residue(Real *g, Real *p, Real *x, int nvar, Real alpha) {
	(void)alpha;
NAMICS_DBG("residue in Newton " << std::endl);
	return std::sqrt(norm2(p,nvar)*norm2(g,nvar)/(1+norm2(x,nvar)));
}

Real SFNewton::linecriterion(Real *g, Real *g0, Real *p, Real *p0, int nvar) {
	(void)p0;
	(void)p;
NAMICS_DBG("linecriterion in Newton " << std::endl);
	Real normg,gg0;
	normg = norm2(g0,nvar);
	gg0 = std::inner_product(g, g + nvar, g0, Real(0));

	gg0=gg0/normg/normg;
	normg = std::pow(norm2(g,nvar)/normg,2);
	if ( (gg0>1 || normg>1) && normg-gg0*std::fabs(gg0)<0.2 ) {
		normg = 1.5*normg;
	}
	if (gg0<0 && normg<2) {
		return 1;
	} else if (normg>10) {
		return .01;
	} else {
		return 0.4*(1+0.75*normg)/(normg-gg0*std::fabs(gg0)+0.1);
	}
}

Real SFNewton::newfunction(Real *g, Real *x, int nvar) {
	(void)x;
NAMICS_DBG("newfunction in Newton " << std::endl);
	return std::pow(norm2(g,nvar),2);
}

void SFNewton::direction(Real *h, Real *p, Real *g, Real *g0, Real *x, int nvar, Real alpha, Real accuracy,bool filter){//done
NAMICS_DBG("direction in Newton " << std::endl);

	newtondirection = true;
	newhessian(h,g,g0,x,p,nvar,accuracy,alpha,filter);
	gausa(h,p,g,nvar);
	gausb(h,p,nvar);
	if (ignore_newton_direction) {
		newtondirection = true;
	} else {
		newtondirection = signdeterminant(h,nvar)>0;
	}
	if ( !newtondirection ) {
		std::transform(p, p + nvar, p, [](Real value) { return -value; });
		if ( e_info && t_info) std::cout << "*";
	}
}

void SFNewton::startderivatives(Real *h, Real *g, Real *x, int nvar){ //done
	(void)x;
NAMICS_DBG("startderivatives in Newton" << std::endl);
	Real diagonal = 1+norm2(g,nvar);
	std::fill_n(h, nvar * nvar, 0);
	for (int i=0; i<nvar; i++) {
		h[i+nvar*i] = diagonal;
	}
}

void SFNewton::resethessian(Real *h,Real *g,Real *x,int nvar){ //done
NAMICS_DBG("resethessian in Newton" << std::endl);
	trouble = 0;
	startderivatives(h,g,x,nvar);
	resetiteration=iterations;
}

void SFNewton::newhessian(Real *h, Real *g, Real *g0, Real *x, Real *p, int nvar, Real accuracy,Real ALPHA,bool filter) {//done
NAMICS_DBG("newhessian in Newton" << std::endl);

	Real dmin=0,sum=0,theta=0,php=0,gg=0,g2=0,py=0,y2=0;
	dmin = 1/std::pow(2.0,nbits); // alternative: DBL_EPSILON or DBL_MIN
	if (!pseudohessian){
		findhessian(h,g,x,nvar,filter);
	} else {
		if (!samehessian && ALPHA!=0 && iterations!=0) {
			std::vector<Real> y(nvar, 0);
			std::vector<Real> hp(nvar, 0);
			py = php = y2 = gg = g2 = 0;
			std::transform(g, g + nvar, g0, y.begin(), std::minus<Real>());
			py = std::inner_product(p, p + nvar, y.begin(), Real(0));
			y2 = std::inner_product(y.begin(), y.end(), y.begin(), Real(0));
			gg = std::inner_product(g, g + nvar, g0, Real(0));
			g2 = std::inner_product(g, g + nvar, g, Real(0));

			if ( !newtondirection ) {
				multiply(hp.data(),1,h,p,nvar);
			} else {
				std::transform(g0, g0 + nvar, hp.begin(), [](Real value) { return -value; });
			}

			php = std::inner_product(p, p + nvar, hp.begin(), Real(0));
			theta = py/(10*dmin+ALPHA*php);

			if ( nvar>=1 && theta>0 && iterations==resetiteration+1 && accuracy > max_accuracy_for_hessian_scaling) {
				if (e_info && hs_info) {
					std::cout << "hessian scaling: " << theta << std::endl;
				}
				ALPHA *= theta;
				py /= theta;
				php /= theta;
				std::transform(p, p + nvar, p, [theta](Real value) { return value / theta; });
				for (int i=0; i<nvar; i++) h[i+nvar*i] *= theta;
			}
			if (nvar>=1) trustfactor *= (4/(std::pow(theta-1,2)+1)+0.5);
			if ( nvar>1 ) {
				sum = ALPHA*std::pow(norm2(p,nvar),2);
				theta = std::fabs(py/(ALPHA*php));
				if ( theta<.01 ) sum /= 0.8;
				else if ( theta>100 ) sum *= theta/50;
				std::transform(y.begin(), y.end(), hp.begin(), y.begin(), [ALPHA](Real yi, Real hpi) { return yi - ALPHA * hpi; });

				updatpos(h,y.data(),p,nvar,1.0/sum);
				trouble -= signdeterminant(h,nvar);
				if ( trouble<0 ) trouble = 0; else if ( trouble>=3 ) resethessian(h,g,x,nvar);
			} else if ( nvar>=1 && py>0 ) {
				trouble = 0;
				theta = py>0.2*ALPHA*php ? 1 : 0.8*ALPHA*php/(ALPHA*php-py);
				if ( theta<1 ) {
					std::transform(y.begin(), y.end(), hp.begin(), y.begin(), [theta, ALPHA](Real yi, Real hpi) {
						return theta * yi + (1 - theta) * ALPHA * hpi;
					});
					py = std::inner_product(p, p + nvar, y.begin(), Real(0));
				}
				updatpos(h,y.data(),y.data(),nvar,1.0/(ALPHA*py));
				updateneg(h,hp.data(),nvar,-1.0/php);
			}
		} else if ( !samehessian ) resethessian(h,g,x,nvar);
	}
}

void SFNewton::numhessian(Real* h,Real* g, Real* x, int nvar,bool filter) {//done
	NAMICS_DBG("numhessian in Newton" << std::endl);
	Real dmax2=0,dmax3=0,di=0;
	std::vector<Real> g1(IV, 0);
	Real xt;
	dmax2 = std::pow(2.0,nbits/2); //alternative 2*std::pow(DBL_EPSILON,-0.5)?
	dmax3 = std::pow(2.0,nbits/3); //alternative 2*std::pow(DBL_EPSILON,-1.0/3)?
	for (int i=0; i<nvar; i++) {
		xt = x[i];
		di = (1/(dmax3*dmax3*std::fabs(h[i+nvar*i])+dmax3+std::fabs(g[i]))
			+1/dmax2)*(1+std::fabs(x[i]));
		x[i] += di;
		COMPUTEG(x,g1.data(),nvar,filter);
		x[i] = xt;
		std::transform(g1.begin(), g1.begin() + nvar, g, &h[nvar * i], [di](Real g1v, Real gv) { return (g1v - gv) / di; });
	}
	COMPUTEG(x,g,nvar,filter);
}

void SFNewton::decomposition(Real *h,int nvar, int &trouble){//done
NAMICS_DBG("decomposition in Newton" << std::endl);
	int ntr=0;
	decompos(h,nvar,ntr);

	if (e_info) {
		if (iterations==0) {
			if (ntr==0) {
				std::cout << "Not too bad.";
			} else {
				if (s_info) std::cout << " sign might be wrong.";
			}
		} else if (ntr>0 && trouble==0) {
			if (ntr==1) {
				if (s_info) std::cout << "Wait a sec.";
			} else {
				if (s_info) std::cout << "some TROUBLES appear.";
			}
		} else if (ntr>trouble+1) {
			for (int i=1; i<= ntr; i++) {
				if (s_info) std::cout << "O";
			}
			if (s_info) std::cout << "H!";
		} else if (trouble>0) {
			if (ntr==0) {
				if (s_info) std::cout << "Here we go.";
			} else if (ntr<trouble-1) {
				if (s_info) std::cout << "Hold on.";
			} else if (ntr < trouble) {
				if (s_info) std::cout << "There is some hope for you.";
			} else if (ntr == trouble) {
				if (s_info) std::cout << "no Progress.";
			} else if (ntr > 4) {
				if (s_info) std::cout << "We won't fix it.";
			} else {
				if (s_info) std::cout << "More troubles.";
			}
		}
		if (iterations==0 || trouble>0 || ntr>0) {
			if (s_info) std::cout <<  std::endl;
		}
	}
	trouble = ntr;

}

void SFNewton::findhessian(Real *h, Real *g, Real *x,int nvar,bool filter) {//done
NAMICS_DBG("findhessian in Newton" << std::endl);
	if ( !samehessian ) {
		if ( iterations==0 ) resethessian(h,g,x,nvar);
		numhessian(h,g,x,nvar,filter); // passes through residuals so check pseudohessian
		if (!pseudohessian) {
			decomposition(h,nvar,trouble);
		}
	}
}


Real SFNewton::newdirection(Real *h, Real *p, Real *p0, Real *g, Real *g0, Real *x, int nvar, Real ALPHA, bool filter) {//done
NAMICS_DBG("newdirection in Newton" << std::endl);

	memcpy(p0, p, sizeof(*p0)*nvar);
	Real accuracy=residue(g,p,x,nvar,ALPHA);
	direction(h,p,g,g0,x,nvar,ALPHA,accuracy,filter);
	return accuracy;
}

void SFNewton::newtrustregion(Real *p0,Real ALPHA_, Real &trustregion, Real& trustfactor, Real delta_max, Real delta_min, int nvar){
NAMICS_DBG("newtrustregion in Newton" << std::endl);
	Real ALPHA=ALPHA_;
	Real normp0 =0;
	normp0= norm2(p0,nvar);

	if ( normp0>0 && trustregion>2*ALPHA*normp0 ) {
		trustregion = 2*ALPHA*normp0;
	}
	trustregion *= trustfactor;
	trustfactor = 1.0;
	if ( trustregion>delta_max ) trustregion = delta_max;
	if ( trustregion<delta_min ) trustregion = delta_min;
}

Real SFNewton::linesearch(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real alphabound,bool filter) {//done
NAMICS_DBG("linesearch in Newton" << std::endl);
	Real newalpha = alphabound<1 ? alphabound : 1;
	newalpha = zero(g,g0,p,x,x0,nvar,newalpha,filter);
	return newalpha;
}

Real SFNewton::zero(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real newalpha,bool filter) {//done
NAMICS_DBG("zero in Newton " << std::endl);
	Real alpha=newalpha;
	bool valid, timedep;
	lineiterations++;
	if (lineiterations==5) {
		memcpy(x, x0, sizeof(*x)*nvar);
		COMPUTEG(x,g,nvar,filter);
		valid = true;
		timedep = true; //to turn off the time-dependence warning which usually is a false one....
		for (int i=0; i<nvar && valid && !timedep; i++) {
			if ( g[i]!=g0[i] && !timedep) {
				std::cout <<"[NEWTON:ERROR?: your functions are time dependent!]"<< std::endl;
				timedep = true;
			} else if (!std::isfinite(g[i]) && valid) {
				std::cout <<"invalid numbers in gradient, reversing search direction"<<std::endl; ;
				valid = false;
				alpha *= -1; // reverse search direction
			}
		}
	}

	std::transform(x0, x0 + nvar, p, x, [alpha](Real x0v, Real pv) { return x0v + alpha * pv; });
	valid = true;

	COMPUTEG(x,g,nvar,filter);
	for (int i=0; i<nvar && valid; i++) {
		if (!std::isfinite(g[i])) {
			valid = false;
			std::cout <<"invalid numbers in gradient"<<std::endl;
				g[i] = 1;
		}
	}
	minimum=newfunction(g,x,nvar);
	return alpha;
}

Real SFNewton::stepchange(Real *g, Real *g0, Real *p, Real *p0, Real *x, Real *x0, int nvar, Real &alpha,bool filter){//done
NAMICS_DBG("stepchange in Newton" << std::endl);
	Real change, crit;
	change = crit = linecriterion(g,g0,p,p0,nvar);
	while ( crit<0.35 && lineiterations<linesearchlimit ) {
		alpha /= 4;
		zero(g,g0,p,x,x0,nvar,alpha,filter);
		crit = linecriterion(g,g0,p,p0,nvar);
		change = 1;
	}
	return change;
}

void SFNewton::COMPUTEG(Real* x, Real* g, int nvar,bool filter) {//done
	int pos=nvar;
	if (filter) {
		for (int i=IV-1; i>=0; i--) {
			if (mask[i]==1) {
				pos--;
				x[i]=x[pos];
			} else x[i]=0;
		}
		residuals(x,g);
		pos=0;
		for (int i=0; i<IV; i++) {
			if (mask[i]==1) {x[pos]=x[i]; g[pos]=g[i];pos++;}
		}
	} else residuals(x,g);
}

void SFNewton::ResetX(Real* x,int nvar,bool filter) { //done
	if (filter) {
		int pos=nvar;
		for (int i=IV-1; i>=0; i--) {
			if (mask[i]==1) {
				pos--;
				x[i]=x[pos];
			} else x[i]=0;
		}
	}
}


bool SFNewton::Message(bool e_info_, bool s_info_, int it_, int iterationlimit_,Real residual_, Real tolerance_, std::string s_) {
	bool e_info=e_info_, s_info=s_info_; std::string s=s_;
	int it=it_, iterationlimit=iterationlimit_;
	Real residual=residual_, tolerance=tolerance_;
	NAMICS_DBG("Message in  Newton " << std::endl);
	bool success=true;
	if (it == iterationlimit) {
		std::cout <<"Warning: "<<s<<"iteration not solved. Residual error= " << residual << std::endl;
		success=false;
	}

  if ((e_info || s_info)) {


		if (e_info) {
			if (it < iterationlimit) std::cout <<s<<"Problem solved." << std::endl;
			if (it < iterationlimit/10) std::cout <<"That was easy." << std::endl;
			if (it > iterationlimit/10 && it < iterationlimit ) std::cout <<"That will do." << std::endl;
			if (it <2 && iterationlimit >1 ) std::cout <<"You hit the nail on the head." << std::endl;
			if (residual > tolerance) { std::cout << " Iterations failed." << std::endl;
				if (residual < tolerance/10) std::cout <<"... I almost made it..." << std::endl;
			}
		}
		if (s_info) std::cout <<it << " iterations used to reach residual " << residual <<"."<< std::endl;
	}
	return success;
}

bool SFNewton::iterate(Real* x,int nvar_,int iterationlimit_,Real tolerance_, Real delta_max_, Real delta_min_,bool filter_) {
NAMICS_DBG("iterate in SFNewton" << std::endl);
	int nvar=nvar_;
	int iterationlimit=iterationlimit_;
	Real tolerance=tolerance_;
	bool success;
	Real delta_max=delta_max_;
	Real delta_min=delta_min_;
	bool filter=filter_;
	std::vector<Real> x0(nvar, 0);
	std::vector<Real> g(nvar, 0);
	std::vector<Real> p(nvar, 0);
	std::vector<Real> p0(nvar, 0);
	std::vector<Real> g0(nvar, 0);
	mask.assign(nvar, 0);
	std::vector<Real> h(nvar * nvar, 0);

	if (nvar<1) {std::cout << "newton has nothing to do; returning the problem" << std::endl; return false;}
	int it=0;
	iterations=it;
	Real alphabound=0;
	Real trustregion=delta_max;
	Real ALPHA=1;
	Real trustfactor =1;

    trouble = resetiteration = 0;
	minAccuracySoFar = 1e30;
	numIterationsSinceHessian=0;

	IV =nvar;
	std::srand (1);
	if (e_info) {std::cout <<"NEWTON has been notified."<< std::endl;
		std::cout << "Your guess:";
	}

	std::copy_n(x, nvar, x0.begin());
	if (filter) {
		std::transform(x, x + nvar, x, [](Real xv) {
			return xv + 1e-10 * (Real)std::rand() / (Real)((unsigned)RAND_MAX + 1);
		});
		residuals(x,g.data());
		std::copy_n(x0.begin(), IV, x);
		int xxx=0;
		for (int i=0; i<nvar; i++) {if (g[i]==0) mask[i]=0; else {xxx++;  mask[i]=1;}}
		nvar=xxx;
		int pos=0;
		for (int i=0; i<IV; i++) {
			if (mask[i]==1) {g[pos]=g[i]; x[pos]=x[i];pos++;}
		}
	}

	newhessian(h.data(),g.data(),g0.data(),x,p.data(),nvar,accuracy,ALPHA,filter);
	minimum=newfunction(g.data(),x,nvar);
	inneriteration(x,g.data(),h.data(),accuracy,delta_max,ALPHA,nvar);
	accuracy=newdirection(h.data(),p.data(),p0.data(),g.data(),g0.data(),x,nvar,ALPHA,filter);
	normg=std::sqrt(minimum);
	accuracy=residue(g.data(),p.data(),x,nvar,ALPHA);

	while ((tolerance < accuracy || tolerance*10<normg) && iterations<iterationlimit && accuracy == std::fabs(accuracy) ) {
		if (e_info)
		if (i_info > 0)
		if (it%i_info == 0) {
#ifdef LongReal
			printf("it =  %i  E = %Le |g| = %Le alpha = %Le \n",it,accuracy,normg,ALPHA);
#else
			printf("it =  %i  E = %e |g| = %e alpha = %e \n",it,accuracy,normg,ALPHA);
#endif
		}
		it++; iterations=it;  lineiterations=0;
		newtrustregion(p0.data(),ALPHA,trustregion,trustfactor,delta_max,delta_min,nvar);  //trustregion and trustfactor are adjusted.
		alphabound = trustregion/(norm2(p.data(),nvar)+1/std::pow(2.0,nbits));
		std::copy_n(x, nvar, x0.begin());
		std::copy_n(g.begin(), nvar, g0.begin());
		ALPHA = linesearch(g.data(),g0.data(),p.data(),x,x0.data(),nvar,alphabound,filter);
		trustfactor *= stepchange(g.data(),g0.data(),p.data(),p0.data(),x,x0.data(),nvar,ALPHA,filter);
		trustfactor *= ALPHA/alphabound;
		inneriteration(x,g.data(),h.data(),accuracy,delta_max,ALPHA,nvar);
		accuracy=newdirection(h.data(),p.data(),p0.data(),g.data(),g0.data(),x,nvar,ALPHA,filter);
		normg=std::sqrt(minimum);
	}
#ifdef LongReal
	if (e_info) printf("it =  %i  E = %Le |g| = %Le alpha = %Le \n",it,accuracy,normg,ALPHA);
#else
	if (e_info) printf("it =  %i  E = %e |g| = %e alpha = %e \n",it,accuracy,normg,ALPHA);
#endif
	success=Message(e_info,s_info,it,iterationlimit,accuracy,tolerance,"");
	ResetX(x,nvar,filter);
	mask.clear();
	return success;
}

void SFNewton::Ax(Real* A, Real* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
NAMICS_DBG("Ax in  SFNewton (own svdcmp) " << std::endl);
	if (N <= 1) {
		X[0] = 1;
		return;
	}

	Eigen::MatrixXd M(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			const Real value = A[i * N + j];
			if (!std::isfinite(value)) throw -2;
			M(i, j) = static_cast<double>(value);
		}
	}

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::VectorXd s = svd.singularValues();
	if (!svd.matrixU().allFinite() || !svd.matrixV().allFinite() || !s.allFinite()) throw -3;
	const Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
	Eigen::VectorXd u_proj = svd.matrixU().transpose() * ones;
	Eigen::VectorXd coeff = Eigen::VectorXd::Zero(N);
	for (int i = 0; i < N; ++i) {
		const double sigma = s(i);
		if (std::abs(sigma) > 1e-14) coeff(i) = u_proj(i) / sigma;
	}
	Eigen::VectorXd result = svd.matrixV() * coeff;
	if (!result.allFinite()) throw -3;
	for (int i = 0; i < N; ++i) X[i] = static_cast<Real>(result(i));
}

void SFNewton::DIIS(Real* x, Real* x_x0, Real* xR, Real* Aij, Real* Apij,Real* Ci, int k, int k_diis, int m, int nvar) {
NAMICS_DBG("DIIS in  SFNewton " << std::endl);
  	int posi;

	if (k_diis>m) {
    	k_diis =m;
		for (int i=1; i<m; i++)
      		for (int j=1; j<m; j++)
		    	Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
	}

	for (int i=0; i<k_diis; i++) {
    	posi = k-k_diis+1+i;
    	if (posi<0) {
      		posi +=m;
		}
	  	Real Dvalue;
		Dvalue = std::inner_product(x_x0 + posi * nvar, x_x0 + (posi + 1) * nvar, x_x0 + k * nvar, Real(0));
		Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dvalue + 1e-9;
		// write to (compressed) matrix Apij
		for (int j=0; j<k_diis; j++)
		    Apij[j+k_diis*i] = Aij[j+m*i];
	}


	Ax(Apij,Ci,k_diis);
	for (int i = 0; i < k_diis; ++i) {
		if (!std::isfinite(Ci[i])) throw -5;
	}

	Real normC=0;

	normC = std::accumulate(Ci, Ci + k_diis, Real(0));
	if (!std::isfinite(normC) || std::abs(normC) <= std::numeric_limits<Real>::epsilon()) throw -6;
	std::transform(Ci, Ci + k_diis, Ci, [normC](Real value) { return value / normC; });
	std::fill_n(x, nvar, 0);
	posi = k-k_diis+1;

	if (posi < 0) posi += m;

	auto accumulate = [&](int coeff_idx, int row) {
		const Real* src = xR + row * nvar;
		const Real coeff = Ci[coeff_idx];
		for (int i = 0; i < nvar; ++i) {
			x[i] += coeff * src[i];
		}
	};
	accumulate(0, posi);
	for (int coeff_idx = 1; coeff_idx < k_diis; ++coeff_idx) {
		int row = k - k_diis + 1 + coeff_idx;
		if (row < 0) row += m;
		accumulate(coeff_idx, row);
	}
	for (int i = 0; i < nvar; ++i) {
		if (!std::isfinite(x[i])) throw -5;
	}
}

Real SFNewton::computeresidual(Real* values, int size) {
  Real residual = 0;
  // Compute residual based on maximum error value
  if (max_g == true) {

	Real* H_array;

		H_array = values;

    auto temp_residual = std::minmax_element(H_array, H_array + size);
    if (std::abs(*temp_residual.first) > std::abs(*temp_residual.second)) {
      residual = std::abs(*temp_residual.first);
    } else {
      residual = std::abs(*temp_residual.second);
    }



  } else {
    // Compute residual based on sum of errors
		residual = std::inner_product(values, values + size, values, Real(0));
		residual = std::sqrt(residual);
  }

  return residual;
}


bool SFNewton::iterate_DIIS(Real*x,int nvar_, int m, int iterationlimit,Real tolerance, Real delta_max, int restart_DIIS) {
	NAMICS_DBG("Iterate_DIIS in SFNewton " << std::endl);
	int nvar=nvar_;
	bool success;
	std::vector<Real> Aij(m * m, 0);
	std::vector<Real> Apij(m * m, 0);
	std::vector<Real> Ci(m, 0);
	std::vector<Real> xR(m * nvar, 0);
	std::vector<Real> x_x0(m * nvar, 0);
	std::vector<Real> x0(nvar, 0);
	std::vector<Real> g(nvar, 0);
	iterations=0;
	int k_diis=0;
	int k=0;
	std::copy_n(x, nvar, x0.begin());
  // mol computephi takes long: moltype = monomer

	try {

		residuals(x,g.data());
		for (int i = 0; i < nvar; ++i) {
			if (!std::isfinite(g[i])) throw -5;
		}

		std::transform(x, x + nvar, g.begin(), x, [delta_max](Real xv, Real gv) { return xv - delta_max * gv; });
		for (int i = 0; i < nvar; ++i) {
			if (!std::isfinite(x[i])) throw -5;
		}
		std::transform(x, x + nvar, x0.begin(), x_x0.begin(), std::minus<Real>());
		std::copy_n(x, nvar, xR.begin());
  		residual = computeresidual(g.data(), nvar);
		if (!std::isfinite(residual)) throw -5;

		if (e_info) printf("DIIS has been notified\n");
#ifdef LongReal
		if (e_info) printf("Your guess = %Le \n",residual);
#else
		if (e_info) printf("Your guess = %e \n",residual);
#endif
		while ( residual > tolerance and iterations < iterationlimit) {
			iterations++;
			if (iterations%restart_DIIS==0) {
				k_diis=0; std::cout<<"!";
			}

			std::copy_n(x, nvar, x0.begin());
			residuals(x,g.data());
			for (int i = 0; i < nvar; ++i) {
				if (!std::isfinite(g[i])) throw -5;
			}
			k=iterations % m; k_diis++; //plek voor laatste opslag
			std::transform(x, x + nvar, g.begin(), x, [delta_max](Real xv, Real gv) { return xv - delta_max * gv; });
			for (int i = 0; i < nvar; ++i) {
				if (!std::isfinite(x[i])) throw -5;
			}
			std::copy_n(x, nvar, xR.begin()+k*nvar);
			std::transform(x, x + nvar, x0.begin(), x_x0.begin() + k * nvar, std::minus<Real>());
			DIIS(x,x_x0.data(),xR.data(),Aij.data(),Apij.data(),Ci.data(),k,k_diis,m,nvar);
    			residual = computeresidual(g.data(), nvar);
			if (!std::isfinite(residual)) throw -5;
			if(e_info && iterations%i_info == 0){
#ifdef LongReal
				printf("iterations = %i g = %Le \n",iterations,residual);
#else
				printf("iterations = %i g = %e \n",iterations,residual);
#endif
			}
		}

		residuals(x,g.data());
		for (int i = 0; i < nvar; ++i) {
			if (!std::isfinite(g[i])) throw -5;
		}
		residual = computeresidual(g.data(), nvar);
		if (!std::isfinite(residual)) throw -5;
		success=Message(e_info,true,iterations,iterationlimit,residual,tolerance,"");

	} catch (int error) {

		if (error == -1)
			std::cerr << "Detected GN not larger than 0." << std::endl;
		if (error == -2)
			std::cerr << "Detected invalid numbers in DIIS matrix in Ax." << std::endl;
		if (error == -3)
			std::cerr << "Detected invalid numbers while solving DIIS coefficients." << std::endl;
		if (error == -4)
			std::cerr << "Detected negative phibulk." << std::endl;
		if (error == -5)
			std::cerr << "Detected invalid numbers in DIIS iterate." << std::endl;
		if (error == -6)
			std::cerr << "Detected invalid DIIS coefficient normalization." << std::endl;
		exit(1);
	}
	return success;
}

bool SFNewton::iterate_RF(Real*x, int nvar_,int iterationlimit,Real tolerance, Real delta_max, std::string s) {
	(void)s;
NAMICS_DBG("Iterate_RF in SFNewton " << std::endl);
	int nvar=nvar_;
	bool success;
	std::vector<Real> x0(nvar, 0);
	std::vector<Real> g(nvar, 0);
	Real a=0, b=0, c=0, fa=0, fb=0, fc=0;
	Real res=100.0;
	int k=0,it=0;

	while ((it<iterationlimit) && (std::abs(res)>tolerance)) {
		if (it>0) {
			it=0;
			std::cout <<"restart regular falsi" << std::endl;
		}
		Real x_start=x[0];
		residuals(x,g.data());
		a=1;
		fa=g[0];
		x[0]=(a+delta_max)*x_start;

		residuals(x,g.data());
		b=x[0]/x_start;

		fb=g[0];
		if(fa==fb) std::cout << "WARNING: The Denominator in Regula Falsi is zero for finding the closest root."<<std::endl;
		c = a - 0.5*((fa*(a-b))/(fa-fb));
		x[0]=c*x_start;
		residuals(x,g.data());
		fc=g[0]; res=fc;
		k=0;

		while((k<iterationlimit/10) && (std::abs(res)>tolerance)){
			c = a-0.5*((fa*(a-b))/(fa-fb)); x[0]=c*x_start;residuals(x,g.data());fc=g[0];
			if(fc*fb<0){
				b=c; x[0]=b*x_start; residuals(x,g.data()); fb=g[0]; res=fb;
			} else {
				a=c; x[0]=a*x_start; residuals(x,g.data()); fa=g[0]; res=fa;
			}
			k++; it++;
			if(fa==fb) std::cout << "WARNING: The Denominator in Regula Falsi is zero for finding the closest root."<<std::endl;
			if (e_info || it>49) {
				if (it==1||it==2||it==4||it==8||it==16||it==32 || (it>50 &&it%10==0)) std::cout << "s_it = " << it << " g = " <<	res << std::endl;
			}
		}
	}

	success=Message(e_info,s_info,it,iterationlimit,residual,tolerance,"");
	return success;
}
