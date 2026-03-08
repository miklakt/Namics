#include <iostream>
#include "LG1Planar.h"

LG1Planar::LG1Planar(const Input& In_,const string& name_): LGrad1(In_,name_) {}

LG1Planar::~LG1Planar() {
NAMICS_DBG("LG1Planar destructor " << endl);}

void LG1Planar:: ComputeLambdas() {
	for (int i=1; i<MX+1; i++) L[i]=1;

	if (fjc>1) {
		for (int i = 0; i < M; i++) {
			L[i] = 1.0/fjc;
			LAMBDA[i] = 1.0/(2*(FJC-1));
			LAMBDA[i + (FJC-1)*M] = 1.0/(2*(FJC-1));
			LAMBDA[i+(FJC-1)/2*M] = 1.0/(FJC-1);
			for (int j = 1; j < FJC/2; j++) {
				LAMBDA[i+j*M] = 1.0/(FJC-1);
				LAMBDA[i+(FJC-j-1)*M] = 1.0/(FJC-1);
			}
		}
	}
}

void LG1Planar::Side(Real *X_side, Real *X, int M) {
NAMICS_DBG(" Side in LG1Planar " << endl);	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0); //set_bounds(X);
	int kk;

	if (fcc_sites) {
		for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i];
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i];
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (1.0/3.0);
	} else {
		if (fjc==1) {
			for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (lambda) * (X)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (lambda) * (X+1)[__i];
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (1.0-2.0*lambda) * (X)[__i];
		} else {
			for (int j = 0; j < FJC/2; j++) {
				kk = (FJC-1)/2-j;
				for (int __i = 0; __i < (M-kk); ++__i) (X_side+kk)[__i] += (X)[__i] * (LAMBDA+j*M+kk)[__i];
				for (int __i = 0; __i < (M-kk); ++__i) (X_side)[__i] += (X+kk)[__i] * (LAMBDA+(FJC-j-1)*M)[__i];
			}
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i] * (LAMBDA+(FJC-1)/2*M)[__i];
		}
	}
}

void LG1Planar::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
NAMICS_DBG(" propagateF in LG1Planar " << endl);
	Real *gs = G+M*FJC*(s_to), *gs_1 = G+M*FJC*(s_from);
	Real *g = G1;

	std::fill_n(gs, M*FJC, 0);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	if (lattice_type==simple_cubic) {
		Real *gz0 = gs_1, *gz1 = gs_1+M, *gz2 = gs_1+2*M;
		Real *gx0 = gs, *gx1 = gs+M, *gx2 = gs+2*M;

		for (int __i = 0; __i < (M-1); ++__i) (gx0+1)[__i] += (P[0]) * (gz0)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gx0+1)[__i] += (4*P[1]) * (gz1)[__i];

		for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (P[1]) * (gz0)[__i];
		for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (2*P[1]+P[0]) * (gz1)[__i];
		for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (P[1]) * (gz2)[__i];

		for (int __i = 0; __i < (M-1); ++__i) (gx2)[__i] += (4*P[1]) * (gz1+1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gx2)[__i] += (P[0]) * (gz2+1)[__i];

		for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
	} else {
		int a,b; Real c;

		for (int p=0; p<FJC; p++){
			a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
			for (int q=0; q<FJC; q++) {
				c=P[abs(-p+q)];
				if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
				if (c!=0) for (int __i = 0; __i < (M-a-b); ++__i) (gs+p*M+a)[__i] += (c) * (gs_1+q*M+b)[__i];
			}
		}
		for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
	}
}

void LG1Planar::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
NAMICS_DBG(" propagateB in LG1Planar " << endl);
	Real *gs = G+M*FJC*(s_to), *gs_1 = G+M*FJC*(s_from);
	Real *g = G1;

	std::fill_n(gs, M*FJC, 0);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	if (lattice_type==simple_cubic) {
		Real *gz0 = gs_1, *gz1 = gs_1+M, *gz2 = gs_1+2*M;
		Real *gx0 = gs,   *gx1 = gs+M,   *gx2 = gs+2*M;

		for (int __i = 0; __i < (M-1); ++__i) (gx0)[__i] += (P[0]) * (gz0+1)[__i];
		for (int __i = 0; __i < (M); ++__i) (gx0)[__i] += (4*P[1]) * (gz1)[__i];

		for (int __i = 0; __i < (M-1); ++__i) (gx1)[__i] += (P[1]) * (gz0+1)[__i];
		for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (2*P[1]+P[0]) * (gz1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gx1+1)[__i] += (P[1]) * (gz2)[__i];

		for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (4*P[1]) * (gz1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gx2+1)[__i] += (P[0]) * (gz2)[__i];

		for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
	} else {
		int a,b; Real c;

		for (int q=FJC-1; q>-1; q--){
			a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
			for (int p=FJC-1; p>-1; p--) {
				c=P[abs(-p+q)];
				if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
				if (c!=0) for (int __i = 0; __i < (M-a-b); ++__i) (gs+p*M+a)[__i] += (c) * (gs_1+q*M+b)[__i];
			}
		}
		for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
	}

}

void LG1Planar::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
NAMICS_DBG(" propagate in LG1Planar " << endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int kk;
	int j;
	std::fill_n(gs, M, 0); set_bounds(gs_1);

	if (fjc==1) {
		for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (lambda) * (gs_1)[__i];
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (1.0-2.0*lambda) * (gs_1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (lambda) * (gs_1+1)[__i];
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
	} else {
		for (j = 0; j < FJC/2; j++) {
			kk = (FJC-1)/2-j;
			for (int __i = 0; __i < (M-kk); ++__i) (gs+kk)[__i] += (gs_1)[__i] * (LAMBDA+j*M+kk)[__i];
			for (int __i = 0; __i < (M-kk); ++__i) (gs)[__i] += (gs_1+kk)[__i] * (LAMBDA+(FJC-j-1)*M)[__i];
		}
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (gs_1)[__i] * (LAMBDA+(FJC-1)/2*M)[__i];
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
	}
}


void LG1Planar::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	std::fill_n(EE, M, 0);
	Real Exmin,Explus;

	pf =pf/2.0*fjc*fjc;
	Explus=psi[fjc-1]-psi[fjc];
	Explus *=Explus;

	for (int x=fjc; x<MX+fjc; x++) {
		Exmin=Explus;
		Explus=psi[x]-psi[x+1];
		Explus *=Explus;
		EE[x]=pf*(Exmin+Explus);
	}
}


void LG1Planar::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, Real* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	Real a,b,c,a_,b_,c_;
	Real epsXplus, epsXmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	C=C*2.0/fjc/fjc;
	epsXplus=eps[fjc-1]+eps[fjc];
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) {
		epsXmin=epsXplus;
		epsXplus=eps[x]+eps[x+1];
		//X[x]=(epsXmin*a  +C*q[x] + epsXplus*c)/(epsXmin+epsXplus);
		if (x==fjc) a=psi[fjc-1]; else a=X[x-1]; //upwind
		X[x]=(epsXmin*a  +C*q[x] + epsXplus*psi[x+1])/(epsXmin+epsXplus);
	}
	for (int __i = 0; __i < (M); ++__i) (g)[__i] = (g)[__i] - (X)[__i];
   } else { //fixedPsi0 is true;
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) {
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 0) psi[x]=0.5*(a+c)+q[x]*C/eps[x];
	}
	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc];a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (int x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 0) {
				psi[x]+=0.25*(c_-a_)*(c-a)/eps[x]*fjc*fjc;
			}
		}
	}
	for (int x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 0) {
		g[x]-=psi[x];
	}
   }
}


void LG1Planar::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, Real* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	Real a,b,c,a_,b_,c_;

	Real C = -e*e/(eps0*k_BT*bond_length);
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) { //for all geometries
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 1) q[x] = -0.5*(a-2*b+c)*fjc*fjc*eps[x]/C;
	}

	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc]; a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (int x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 1) q[x]-=0.25*(c_-a_)*(c-a)*fjc*fjc/C;
		}
	}
	for (int x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 1) {
		g[x]=-q[x];
	}
}
bool LG1Planar:: PutMask(Real* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	(void)R;
	(void)pz;
	(void)py;
	(void)px;
	(void)MASK;
	bool success=true;
	cout <<"PutMask does not make sence in planar 1 gradient system " << endl;
	return success;
}

Real LG1Planar::MomentPlanar(Real* X, int n, Real Z0){
	Real result=0;
	for (int z=0; z<M; z++) {
		result +=X[z]*pow(z-Z0,n);
	}
	return result;
}

