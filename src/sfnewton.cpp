#include <iostream>
#include "sfnewton.h"
#include "tools_host.h"

using namespace std;

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
Addaptation for using vector class in namics (filtering with mask):
Frans Leermakers, Wageningen Agricultural University, NL.

C Copyright (2018) Wageningen University, NL.

 *NO PART OF THIS WORK MAY BE REPRODUCED, EITHER ELECTRONICALLY OF OTHERWISE*

*/
       nbits = numeric_limits<Real>::digits;	//nbits=52;
       linesearchlimit = iterations=lineiterations=numIterationsSinceHessian = trouble=resetiteration=0;
       trouble = resetiteration = 0;
	numReverseDirection =0;
	trustregion =0.0;
	pseudohessian = samehessian = false;
	d_info = e_info = g_info = h_info = s_info = x_info = false;
	newtondirection  = false ;
	ignore_newton_direction = true;
	i_info=1;
	max_accuracy_for_hessian_scaling = 0.1;
	linesearchlimit = 20;
	linetolerance = 9e-1;
	epsilon = 0.1/pow(2.0,nbits/2);
	minAccuracySoFar = 1e30;
	reverseDirectionRange = 50;
	resetHessianCriterion = 1e5;
	reset_pseudohessian = false;
	accuracy=1e30;
	numIterationsSinceHessian=0;
	maxFrReverseDirection=0.4;
	numIterationsForHessian=100;
	minAccuracyForHessian=0.1;
	reverseDirection=(int*) malloc(reverseDirectionRange*sizeof(int)); H_Zero(reverseDirection,reverseDirectionRange);

}

SFNewton::~SFNewton() {
  free(reverseDirection);
}


void SFNewton::multiply(Real *v,Real alpha, Real *h, Real *w, int nvar) { //done
if(debug) cout <<"multiply in Newton" << endl;
	int i=0,i1=0,j=0;
	Real sum=0;
	Real *x = new Real[nvar];
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
	delete [] x;
}


Real SFNewton::norm2(Real*x, int nvar) { //done
if(debug) cout <<"norm2 in Newton" << endl;
	Real sum=0;
	for (int i=0; i<nvar; i++) sum += pow(x[i],2);
	return sqrt(sum);
}

int SFNewton::signdeterminant(Real*h,int nvar) { //dome
if(debug) cout <<"signdeterminant in Newton" << endl;
	int sign=1;
	for (int i=0; i<nvar; i++) {
		if ( h[i+i*nvar]<0 ) {
			sign = -sign;
		}
	}
	return sign;
}

void SFNewton::updateneg(Real *l,Real *w, int nvar, Real alpha) { //done
if(debug) cout <<"updateneg in Newton" << endl;
	int i=0,i1=0,j=0;
	Real dmin=0,sum=0,b=0,d=0,p=0,lji=0,t=0;
	dmin = 1.0/pow(2.0,54);
	alpha = sqrt(-alpha);
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
if(debug) cout <<"decompos in Newton" << endl;
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
if(debug) cout <<"updatepos in Newton" << endl;
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
if(debug) cout <<"gausa in Newton" << endl;
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
if(debug) cout <<"gausb in Newton " << endl;
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
if(debug) cout <<"residue in Newton " << endl;
	return sqrt(norm2(p,nvar)*norm2(g,nvar)/(1+norm2(x,nvar)));
}

Real SFNewton::linecriterion(Real *g, Real *g0, Real *p, Real *p0, int nvar) {
	(void)p0;
	(void)p;
if(debug) cout <<"linecriterion in Newton " << endl;
	Real normg,gg0;
	normg = norm2(g0,nvar);
  	(gg0) = 0; for (int __i = 0; __i < (nvar); ++__i) (gg0) += (g)[__i] * (g0)[__i];

	gg0=gg0/normg/normg;
	normg = pow(norm2(g,nvar)/normg,2);
	if ( (gg0>1 || normg>1) && normg-gg0*fabs(gg0)<0.2 ) {
		normg = 1.5*normg;
	}
	if (gg0<0 && normg<2) {
		return 1;
	} else if (normg>10) {
		return .01;
	} else {
		return 0.4*(1+0.75*normg)/(normg-gg0*fabs(gg0)+0.1);
	}
}

Real SFNewton::newfunction(Real *g, Real *x, int nvar) {
	(void)x;
if(debug) cout <<"newfunction in Newton " << endl;
	return pow(norm2(g,nvar),2);
}

void SFNewton::direction(Real *h, Real *p, Real *g, Real *g0, Real *x, int nvar, Real alpha, Real accuracy,bool filter){//done
if(debug) cout <<"direction in Newton " << endl;

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
		for (int i=0; i<nvar; i++) {
			p[i] *= -1;
		}
		if ( e_info && t_info) cout << "*";
	}
}

void SFNewton::startderivatives(Real *h, Real *g, Real *x, int nvar){ //done
	(void)x;
if(debug) cout <<"startderivatives in Newton" << endl;
	Real diagonal = 1+norm2(g,nvar);
	H_Zero(h,nvar*nvar);
	for (int i=0; i<nvar; i++) {
		h[i+nvar*i] = diagonal;
	}
}

void SFNewton::resethessian(Real *h,Real *g,Real *x,int nvar){ //done
if(debug) cout <<"resethessian in Newton" << endl;
	trouble = 0;
	startderivatives(h,g,x,nvar);
	resetiteration=iterations;
}

void SFNewton::newhessian(Real *h, Real *g, Real *g0, Real *x, Real *p, int nvar, Real accuracy,Real ALPHA,bool filter) {//done
if(debug) cout <<"newhessian in Newton" << endl;

	Real dmin=0,sum=0,theta=0,php=0,dg=0,gg=0,g2=0,py=0,y2=0;
	dmin = 1/pow(2.0,nbits); // alternative: DBL_EPSILON or DBL_MIN
	if (!pseudohessian){
		findhessian(h,g,x,nvar,filter);
	} else {
		if (!samehessian && ALPHA!=0 && iterations!=0) {
			Real *y = new Real[nvar];
			Real *hp = new Real[nvar];
			py = php = y2 = gg = g2 = 0;
			for (int i=0; i<nvar; i++) {
				y[i] = dg = g[i]-g0[i];
				py += p[i]*dg;
				y2 += pow(y[i],2);
				gg += g[i]*g0[i];
				g2 += pow(g[i],2);
			}

			if ( !newtondirection ) {
				multiply(hp,1,h,p,nvar);
			} else {
				for (int i=0; i<nvar; i++) hp[i] = -g0[i];
			}

			(php) = 0; for (int __i = 0; __i < (nvar); ++__i) (php) += (p)[__i] * (hp)[__i];
			theta = py/(10*dmin+ALPHA*php);

			if ( nvar>=1 && theta>0 && iterations==resetiteration+1 && accuracy > max_accuracy_for_hessian_scaling) {
				if (e_info && hs_info) {
					cout << "hessian scaling: " << theta << endl;
				}
				ALPHA *= theta;
				py /= theta;
				php /= theta;
				for (int i=0; i<nvar; i++) {
					p[i] /= theta;
					h[i+nvar*i] *= theta;
				}
			}
			if (nvar>=1) trustfactor *= (4/(pow(theta-1,2)+1)+0.5);
			if ( nvar>1 ) {
				sum = ALPHA*pow(norm2(p,nvar),2);
				theta = fabs(py/(ALPHA*php));
				if ( theta<.01 ) sum /= 0.8;
				else if ( theta>100 ) sum *= theta/50;
				for (int i=0; i<nvar; i++) y[i] -= ALPHA*hp[i];

				updatpos(h,y,p,nvar,1.0/sum);
				trouble -= signdeterminant(h,nvar);
				if ( trouble<0 ) trouble = 0; else if ( trouble>=3 ) resethessian(h,g,x,nvar);
			} else if ( nvar>=1 && py>0 ) {
				trouble = 0;
				theta = py>0.2*ALPHA*php ? 1 : 0.8*ALPHA*php/(ALPHA*php-py);
				if ( theta<1 ) {
					py = 0;
					for (int i=0; i<nvar; i++) {
						y[i] = theta*y[i]+(1-theta)*ALPHA*hp[i];
						py += p[i]*y[i];
					}
				}
				updatpos(h,y,y,nvar,1.0/(ALPHA*py));
				updateneg(h,hp,nvar,-1.0/php);
			}

			delete [] y;
			delete [] hp;
		} else if ( !samehessian ) resethessian(h,g,x,nvar);
	}
}

void SFNewton::numhessian(Real* h,Real* g, Real* x, int nvar,bool filter) {//done
if(debug) cout <<"numhessian in Newton" << endl;
	Real dmax2=0,dmax3=0,di=0;
	Real *g1;
	g1 = new Real[IV];
	Real xt;
	dmax2 = pow(2.0,nbits/2); //alternative 2*pow(DBL_EPSILON,-0.5)?
	dmax3 = pow(2.0,nbits/3); //alternative 2*pow(DBL_EPSILON,-1.0/3)?
	for (int i=0; i<nvar; i++) {
		xt = x[i];
		di = (1/(dmax3*dmax3*fabs(h[i+nvar*i])+dmax3+fabs(g[i]))
			+1/dmax2)*(1+fabs(x[i]));
		//}
		x[i] += di;
		COMPUTEG(x,g1,nvar,filter);
		x[i] = xt;
		for (int j=0; j<nvar; j++ ) {
			h[j+nvar*i] = (g1[j]-g[j])/di;
		}
	}
	delete [] g1;
	COMPUTEG(x,g,nvar,filter);
}

void SFNewton::decomposition(Real *h,int nvar, int &trouble){//done
if(debug) cout <<"decomposition in Newton" << endl;
	int ntr=0;
	decompos(h,nvar,ntr);

	if (e_info) {
		if (iterations==0) {
			if (ntr==0) {
				cout << "Not too bad.";
			} else {
				if (s_info) cout << " sign might be wrong.";
			}
		} else if (ntr>0 && trouble==0) {
			if (ntr==1) {
				if (s_info) cout << "Wait a sec.";
			} else {
				if (s_info) cout << "some TROUBLES appear.";
			}
		} else if (ntr>trouble+1) {
			for (int i=1; i<= ntr; i++) {
				if (s_info) cout << "O";
			}
			if (s_info) cout << "H!";
		} else if (trouble>0) {
			if (ntr==0) {
				if (s_info) cout << "Here we go.";
			} else if (ntr<trouble-1) {
				if (s_info) cout << "Hold on.";
			} else if (ntr < trouble) {
				if (s_info) cout << "There is some hope for you.";
			} else if (ntr == trouble) {
				if (s_info) cout << "no Progress.";
			} else if (ntr > 4) {
				if (s_info) cout << "We won't fix it.";
			} else {
				if (s_info) cout << "More troubles.";
			}
		}
		if (iterations==0 || trouble>0 || ntr>0) {
			if (s_info) cout <<  endl;
		}
	}
	trouble = ntr;

}

void SFNewton::findhessian(Real *h, Real *g, Real *x,int nvar,bool filter) {//done
if(debug) cout <<"findhessian in Newton" << endl;
	if ( !samehessian ) {
		if ( iterations==0 ) resethessian(h,g,x,nvar);
		numhessian(h,g,x,nvar,filter); // passes through residuals so check pseudohessian
		if (!pseudohessian) {
			decomposition(h,nvar,trouble);
		}
	}
}


Real SFNewton::newdirection(Real *h, Real *p, Real *p0, Real *g, Real *g0, Real *x, int nvar, Real ALPHA, bool filter) {//done
if(debug) cout <<"newdirection in Newton" << endl;

	memcpy(p0, p, sizeof(*p0)*nvar);
	Real accuracy=residue(g,p,x,nvar,ALPHA);
	direction(h,p,g,g0,x,nvar,ALPHA,accuracy,filter);
	return accuracy;
}

void SFNewton::newtrustregion(Real *p0,Real ALPHA_, Real &trustregion, Real& trustfactor, Real delta_max, Real delta_min, int nvar){
if(debug) cout <<"newtrustregion in Newton" << endl;
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
if(debug) cout <<"linesearch in Newton" << endl;
	Real newalpha = alphabound<1 ? alphabound : 1;
	newalpha = zero(g,g0,p,x,x0,nvar,newalpha,filter);
	return newalpha;
}

Real SFNewton::zero(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real newalpha,bool filter) {//done
if(debug) cout <<"zero in Newton " << endl;
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
				cout <<"[NEWTON:ERROR?: your functions are time dependent!]"<< endl;
				timedep = true;
			} else if (!isfinite(g[i]) && valid) {
				cout <<"invalid numbers in gradient, reversing search direction"<<endl; ;
				valid = false;
				alpha *= -1; // reverse search direction
			}
		}
	}

	for (int i=0; i<nvar; i++) x[i] = x0[i]+alpha*p[i];
	valid = true;

	COMPUTEG(x,g,nvar,filter);
	for (int i=0; i<nvar && valid; i++) {
		if (!isfinite(g[i])) {
			valid = false;
			cout <<"invalid numbers in gradient"<<endl;
				g[i] = 1;
		}
	}
	minimum=newfunction(g,x,nvar);
	return alpha;
}

Real SFNewton::stepchange(Real *g, Real *g0, Real *p, Real *p0, Real *x, Real *x0, int nvar, Real &alpha,bool filter){//done
if(debug) cout <<"stepchange in Newton" << endl;
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


bool SFNewton::Message(bool e_info_, bool s_info_, int it_, int iterationlimit_,Real residual_, Real tolerance_, string s_) {
	bool e_info=e_info_, s_info=s_info_; string s=s_;
	int it=it_, iterationlimit=iterationlimit_;
	Real residual=residual_, tolerance=tolerance_;
	if (debug) cout <<"Message in  Newton " << endl;
	bool success=true;
	if (it == iterationlimit) {
		cout <<"Warning: "<<s<<"iteration not solved. Residual error= " << residual << endl;
		success=false;
	}

  if ((e_info || s_info)) {


		if (e_info) {
			if (it < iterationlimit) cout <<s<<"Problem solved." << endl;
			if (it < iterationlimit/10) cout <<"That was easy." << endl;
			if (it > iterationlimit/10 && it < iterationlimit ) cout <<"That will do." << endl;
			if (it <2 && iterationlimit >1 ) cout <<"You hit the nail on the head." << endl;
			if (residual > tolerance) { cout << " Iterations failed." << endl;
				if (residual < tolerance/10) cout <<"... I almost made it..." << endl;
			}
		}
		if (s_info) cout <<it << " iterations used to reach residual " << residual <<"."<< endl;
	}
	return success;
}

bool SFNewton::iterate(Real* x,int nvar_,int iterationlimit_,Real tolerance_, Real delta_max_, Real delta_min_,bool filter_) {
if(debug) cout <<"iterate in SFNewton" << endl;
	int nvar=nvar_;
	int iterationlimit=iterationlimit_;
	Real tolerance=tolerance_;
	bool success;
	Real delta_max=delta_max_;
	Real delta_min=delta_min_;
	bool filter=filter_;
  Real* x0 = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(x0, nvar, 0);
  Real* g = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(g, nvar, 0);
  Real* p = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(p, nvar, 0);
  Real* p0 = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(p0, nvar, 0);
  Real* g0  = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(g0, nvar, 0);
  mask = (int*) malloc(nvar*sizeof(int));
  Real* h = (Real*) malloc(nvar*nvar*sizeof(Real)); std::fill_n(h, nvar*nvar, 0);

	if (nvar<1) {cout << "newton has nothing to do; returning the problem" << endl; return false;}
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
	srand (1);
	if (e_info) {cout <<"NEWTON has been notified."<< endl;
		cout << "Your guess:";
	}

	std::copy_n(x, nvar, x0);
	if (filter) {
		for (int i=0; i<nvar; i++) x[i]+=1e-10*(Real)rand() / (Real)((unsigned)RAND_MAX + 1);
		residuals(x,g);
		std::copy_n(x0, IV, x);
		int xxx=0;
		for (int i=0; i<nvar; i++) {if (g[i]==0) mask[i]=0; else {xxx++;  mask[i]=1;}}
		nvar=xxx;
		int pos=0;
		for (int i=0; i<IV; i++) {
			if (mask[i]==1) {g[pos]=g[i]; x[pos]=x[i];pos++;}
		}
	}

	newhessian(h,g,g0,x,p,nvar,accuracy,ALPHA,filter);
	minimum=newfunction(g,x,nvar);
	inneriteration(x,g,h,accuracy,delta_max,ALPHA,nvar);
	accuracy=newdirection(h,p,p0,g,g0,x,nvar,ALPHA,filter);
	normg=sqrt(minimum);
	accuracy=residue(g,p,x,nvar,ALPHA);

	while ((tolerance < accuracy || tolerance*10<normg) && iterations<iterationlimit && accuracy == fabs(accuracy) ) {
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
		newtrustregion(p0,ALPHA,trustregion,trustfactor,delta_max,delta_min,nvar);  //trustregion and trustfactor are adjusted.
		alphabound = trustregion/(norm2(p,nvar)+1/pow(2.0,nbits));
		std::copy_n(x, nvar, x0);
		std::copy_n(g, nvar, g0);
		ALPHA = linesearch(g,g0,p,x,x0,nvar,alphabound,filter);
		trustfactor *= stepchange(g,g0,p,p0,x,x0,nvar,ALPHA,filter);
		trustfactor *= ALPHA/alphabound;
		inneriteration(x,g,h,accuracy,delta_max,ALPHA,nvar);
		accuracy=newdirection(h,p,p0,g,g0,x,nvar,ALPHA,filter);
		normg=sqrt(minimum);
	}
#ifdef LongReal
	if (e_info) printf("it =  %i  E = %Le |g| = %Le alpha = %Le \n",it,accuracy,normg,ALPHA);
#else
	if (e_info) printf("it =  %i  E = %e |g| = %e alpha = %e \n",it,accuracy,normg,ALPHA);
#endif
	success=Message(e_info,s_info,it,iterationlimit,accuracy,tolerance,"");
	ResetX(x,nvar,filter);
 	free(x0);
	free(g);
	free(p);
	free(p0);
	free(g0);
	free(h);
	free(mask);

  mask = NULL;
	return success;
}


bool SFNewton::iterate_Picard(Real* x,int nvar, int iterationlimit, Real tolerance, Real delta_max) {
if(debug) cout <<"Iterate_Picard in  SFNewton " << endl;

Real* h  = (Real*) malloc(sizeof(Real));
Real* g = (Real*) malloc(nvar*sizeof(Real));

	bool success=true;
	int it;
	Real residual;

	it=0;
	if (e_info) {cout <<"Picard has been notified" << endl;
		cout << "Your guess:";
	}
	residuals(x,g);
	residual=computeresidual(g,nvar);
	while (residual > tolerance && it < iterationlimit) {
		if(it%i_info == 0){
#ifdef LongReal
			printf("it = %i g = %Le \n",it,residual);
#else
			printf("it = %i g = %e \n",it,residual);
#endif
		}
		for (int __i = 0; __i < (nvar); ++__i) (x)[__i] += (delta_max) * (g)[__i];
		residual=computeresidual(g,nvar);
		//inneriteration(x,g,h,residual,nvar);
		it++;
	}
	success=Message(e_info,s_info,it,iterationlimit,residual,tolerance,"");
free(h); free(g);
	return success;
}

void SFNewton::Ax(Real* A, Real* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
if(debug) cout <<"Ax in  SFNewton (own svdcmp) " << endl;
	if (N <= 1) {
		X[0] = 1;
		return;
	}

	// Solve A * X = 1 using Gaussian elimination with partial pivoting.
	vector<Real> M(A, A + N * N);
	vector<Real> b(N, 1.0);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (M[i * N + j] != M[i * N + j]) throw -2;
		}
	}

	for (int col = 0; col < N; ++col) {
		int pivot = col;
		Real max_abs = fabs(M[col * N + col]);
		for (int row = col + 1; row < N; ++row) {
			const Real v = fabs(M[row * N + col]);
			if (v > max_abs) {
				max_abs = v;
				pivot = row;
			}
		}
		if (max_abs < 1e-14) throw -3;

		if (pivot != col) {
			for (int j = col; j < N; ++j) {
				std::swap(M[col * N + j], M[pivot * N + j]);
			}
			std::swap(b[col], b[pivot]);
		}

		const Real diag = M[col * N + col];
		for (int row = col + 1; row < N; ++row) {
			const Real factor = M[row * N + col] / diag;
			if (factor == 0) continue;
			M[row * N + col] = 0;
			for (int j = col + 1; j < N; ++j) {
				M[row * N + j] -= factor * M[col * N + j];
			}
			b[row] -= factor * b[col];
		}
	}

	for (int row = N - 1; row >= 0; --row) {
		Real sum = b[row];
		for (int j = row + 1; j < N; ++j) {
			sum -= M[row * N + j] * X[j];
		}
		const Real diag = M[row * N + row];
		if (fabs(diag) < 1e-14) throw -3;
		X[row] = sum / diag;
	}
}

void SFNewton::DIIS(Real* x, Real* x_x0, Real* xR, Real* Aij, Real* Apij,Real* Ci, int k, int k_diis, int m, int nvar) {
if(debug) cout <<"DIIS in  SFNewton " << endl;
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
    		(Dvalue) = 0; for (int __i = 0; __i < (nvar); ++__i) (Dvalue) += (x_x0+posi*nvar)[__i] * (x_x0+k*nvar)[__i];
		Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dvalue + 1e-9;
		// write to (compressed) matrix Apij
		for (int j=0; j<k_diis; j++)
		    Apij[j+k_diis*i] = Aij[j+m*i];
	}


	Ax(Apij,Ci,k_diis);

	Real normC=0;

	for (int i=0; i<k_diis; i++)
    	normC +=Ci[i];
	for (int i=0; i<k_diis; i++)
    	Ci[i] =Ci[i]/normC;
	std::fill_n(x, nvar, 0);
	posi = k-k_diis+1;

  	if (posi<0)
    	posi +=m;

	Xr_times_ci(posi, k_diis, k, m, nvar, x, xR, d_Ci);
}

Real SFNewton::computeresidual(Real* array, int size) {
  Real residual = 0;
  // Compute residual based on maximum error value
  if (max_g == true) {

	Real* H_array;

		H_array = array;

    auto temp_residual = minmax_element(H_array, H_array+size);
    if(abs(*temp_residual.first) > abs(*temp_residual.second) ) {
      residual = abs(*temp_residual.first);
    } else {
      residual = abs(*temp_residual.second);
    }



  } else {
    // Compute residual based on sum of errors
		residual = H_Dot(array,array,size);
		residual = sqrt(residual);
  }

  return residual;
}


bool SFNewton::iterate_BRR(Real*x,int nvar_, int m, int iterationlimit,Real tolerance, Real delta_max) {
	(void)x;
	(void)nvar_;
	(void)m;
	(void)iterationlimit;
	(void)tolerance;
	(void)delta_max;
	cout << "BRR method is disabled in this minimal build. Use DIIS or pseudohessian." << endl;
	return false;
}

bool SFNewton::iterate_DIIS(Real*x,int nvar_, int m, int iterationlimit,Real tolerance, Real delta_max, int restart_DIIS) {
if(debug) cout <<"Iterate_DIIS in SFNewton " << endl;
	int nvar=nvar_;
	bool success;
  Real* Aij = (Real*) malloc(m*m*sizeof(Real)); H_Zero(Aij,m*m);
  Real* Apij = (Real*) malloc(m*m*sizeof(Real)); H_Zero(Apij,m*m);
  Real* Ci = (Real*) malloc(m*sizeof(Real)); H_Zero(Ci,m);
  d_Ci = Ci;
  Real* xR = (Real*) malloc(m*nvar*sizeof(Real)); std::fill_n(xR, m*nvar, 0);
  Real* x_x0 = (Real*) malloc(m*nvar*sizeof(Real)); std::fill_n(x_x0, m*nvar, 0);
  Real* x0 = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(x0, nvar, 0);
  Real* g = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(g, nvar, 0);
	iterations=0;
	int k_diis=0;
	int k=0;
	std::copy_n(x, nvar, x0);
  // mol computephi takes long: moltype = monomer

	try {

		residuals(x,g);

		for (int __i = 0; __i < (nvar); ++__i) (x)[__i] += (-delta_max) * (g)[__i];
		for (int __i = 0; __i < (nvar); ++__i) (x_x0)[__i] = (x)[__i] - (x0)[__i];
		std::copy_n(x, nvar, xR);
  		residual = computeresidual(g, nvar);

		if (e_info) printf("DIIS has been notified\n");
#ifdef LongReal
		if (e_info) printf("Your guess = %Le \n",residual);
#else
		if (e_info) printf("Your guess = %e \n",residual);
#endif
		while ( residual > tolerance and iterations < iterationlimit) {
			iterations++;
			if (iterations%restart_DIIS==0) {
				k_diis=0; cout<<"!";
			}

			std::copy_n(x, nvar, x0);
			residuals(x,g);
			k=iterations % m; k_diis++; //plek voor laatste opslag
			for (int __i = 0; __i < (nvar); ++__i) (x)[__i] += (-delta_max) * (g)[__i];
			std::copy_n(x, nvar, xR+k*nvar);
			for (int __i = 0; __i < (nvar); ++__i) (x_x0+k*nvar)[__i] = (x)[__i] - (x0)[__i];
			DIIS(x,x_x0,xR,Aij,Apij,Ci,k,k_diis,m,nvar);
    			residual = computeresidual(g, nvar);
			if(e_info && iterations%i_info == 0){
#ifdef LongReal
				printf("iterations = %i g = %Le \n",iterations,residual);
#else
				printf("iterations = %i g = %e \n",iterations,residual);
#endif
			}
		}

		success=Message(e_info,true,iterations,iterationlimit,residual,tolerance,"");

	} catch (int error) {

		if (error == -1)
			cerr << "Detected GN not larger than 0." << endl;
		if (error == -2)
			cerr << "Detected nan in U in Ax." << endl;
		if (error == -3)
			cerr << "Detected nan in svdcmp." << endl;
		if (error == -4)
			cerr << "Detected negative phibulk." << endl;
		free(Aij);free(Ci);free(Apij);
  		free(xR);free(x_x0);free(x0);free(g);

		throw error;
	}
  	free(Aij);
	free(Apij);
	free(Ci);
  	free(xR);
	free(x_x0);
	free(x0);
	free(g);
	return success;
}

bool SFNewton::iterate_RF(Real*x, int nvar_,int iterationlimit,Real tolerance, Real delta_max, string s) {
	(void)s;
if(debug) cout <<"Iterate_RF in SFNewton " << endl;
	int nvar=nvar_;
	bool success;
	Real* x0 = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(x0, nvar, 0);
	Real* g = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(g, nvar, 0);
	Real a=0, b=0, c=0, fa=0, fb=0, fc=0;
	Real res=100.0;
	int k=0,it=0;

	while ((it<iterationlimit) && (abs(res)>tolerance)) {
		if (it>0) {
			it=0;
			cout <<"restart regular falsi" << endl;
		}
		Real x_start=x[0];
		residuals(x,g);
		a=1;
		fa=g[0];
		x[0]=(a+delta_max)*x_start;

		residuals(x,g);
		b=x[0]/x_start;

		fb=g[0];
		if(fa==fb) cout << "WARNING: The Denominator in Regula Falsi is zero for finding the closest root."<<endl;
		c = a - 0.5*((fa*(a-b))/(fa-fb));
		x[0]=c*x_start;
		residuals(x,g);
		fc=g[0]; res=fc;
		//it+=k;
		k=0;

		while((k<iterationlimit/10) && (abs(res)>tolerance)){
			c = a-0.5*((fa*(a-b))/(fa-fb)); x[0]=c*x_start;residuals(x,g);fc=g[0];
			if(fc*fb<0){
				b=c; x[0]=b*x_start; residuals(x,g); fb=g[0]; res=fb;
			} else {
				a=c; x[0]=a*x_start; residuals(x,g); fa=g[0]; res=fa;
			}
			k++; it++;
			if(fa==fb) cout << "WARNING: The Denominator in Regula Falsi is zero for finding the closest root."<<endl;
			if (e_info || it>49) {
				if (it==1||it==2||it==4||it==8||it==16||it==32 || (it>50 &&it%10==0)) cout << "s_it = " << it << " g = " <<	res << endl;
			}
		}
	}

	success=Message(e_info,s_info,it,iterationlimit,residual,tolerance,"");
	free(x0);free(g);
	return success;
}



bool SFNewton::iterate_conjugate_gradient(Real *x, int nvar,int iterationlimit , Real tolerance, Real deltamax) {
// Based on An Introduction to the Conjugate Gradient Method Without the Agonizing Pain Edition 1 1/4 - Jonathan Richard Shewchuk
// CG with Newton-Raphson and Fletcher-Reeves
	int  k=0, j=0, j_max =linesearchlimit;
	bool success=true;

	Real inner_err=0, delta_new=0, delta_old=0, delta_d=0, delta_mid = 0;

	Real teller=0, noemer=0;
  Real alpha=0, beta=0;
	Real rd=0;
	bool proceed;
  int iterations=0;
  Real* g = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(g, nvar, 0);
  Real* dg = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(dg, nvar, 0);
  Real* r = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(r, nvar, 0);
  Real* r_old = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(r_old, nvar, 0);
  Real* d = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(d, nvar, 0);
  Real* x0 = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(x0, nvar, 0);
  Real* H_d = (Real*) malloc(nvar*sizeof(Real)); std::fill_n(H_d, nvar, 0);

	if ( e_info ) {
		cout<<"Nonlinear conjugate gradients with Newton-Raphson and Fletcher-Reeves has been notified."<<endl;

	}
	residuals(x,g);

	for (int z=0; z<nvar; z++) {
		r[z] = -g[z];
		d[z]=r[z];
		delta_new += r[z]*r[z];
	}

	accuracy= pow(delta_new,0.5);

	cout << "i = " << iterations << " |g| = "<< accuracy << endl;
	while (tolerance < accuracy && iterations<iterationlimit) {
		j=0;
		delta_d=0;
		for (int z=0; z<nvar; z++) delta_d += d[z]*d[z];
		inner_err = pow(delta_d,0.5);
		proceed=true;
    // line search
		while (proceed) {
			teller=0;

			for (int z=0; z<nvar; z++) {
				teller -= g[z]*d[z];

				x0[z]=x[z];
			}
			Hd(H_d,d,x,x0,g,dg,nvar);
			noemer=0;

			for (int z=0; z<nvar; z++) noemer +=H_d[z]*d[z];
			alpha=teller/noemer;

			for (int z=0; z<nvar; z++) {x[z]=x0[z]+deltamax*alpha*d[z];}
			residuals(x,g);
			j++;
			proceed =(j<j_max && alpha*inner_err> 1e-8);
		}


   		for (int z=0; z<nvar; z++) r_old[z]=r[z];
		for (int z=0; z<nvar; z++) r[z]=-g[z];
		delta_old=delta_new;
		delta_new=0;
    		for (int z=0; z<nvar; z++) delta_mid += r[z]*r_old[z];
		for (int z=0; z<nvar; z++) delta_new += r[z]*r[z];
		beta =  (delta_new-delta_mid)/delta_old;
		for (int z=0; z<nvar; z++) d[z]=r[z]+beta*d[z];
		k++;
		rd=0;
		for (int z=0; z<nvar; z++) rd +=r[z]*d[z];

		if (k == nvar || rd<=0) {
			k=0; beta=0;

			for (int z=0; z<nvar; z++) d[z]=r[z];
		}
		iterations++;
		accuracy = pow(delta_new,0.5);
		if ( e_info ) {
			if (iterations%i_info == 0)
			cout << "i = " << iterations << " |g| = "<< accuracy << "  alpha("<<j<<") = " << alpha << "  beta = " << beta << endl;
		}
	}
	if (iterations==iterationlimit) success=false;
	Message(e_info,true,iterations,iterationlimit,accuracy,tolerance,"");

  free(H_d); free(g);free(dg);free(r);free(r_old);free(d);free(x0);
  return success;
}


void SFNewton::Hd(Real *H_q, Real *q, Real *x, Real *x0, Real *g, Real* dg, Real nvar) {


	Real epsilon = 2e-8; //double precision. Machine error =10^-16; epsilon = 2 sqrt(precision)

	Real normq = norm2(q,nvar);
	Real normx = norm2(x0,nvar);
	Real delta = epsilon* (1+normx)/normq;



	for (int i=0; i<nvar; i++)
    x[i] = x0[i] + delta*q[i];

	residuals(x,dg);
  
	for (int i=0; i<nvar; i++) H_q[i] = (dg[i]-g[i])/delta;
}
