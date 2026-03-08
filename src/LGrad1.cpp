#include <iostream>
#include <string>
#include "lattice.h"
#include "LGrad1.h"

//planar geometry is in LG1Planar.cpp

LGrad1::LGrad1(const Input& In_,const string& name_): Lattice(In_,name_) {
NAMICS_DBG("LGrad1 constructor " << endl);}

LGrad1::~LGrad1() {
NAMICS_DBG("LGrad1 destructor " << endl);}

void LGrad1:: ComputeLambdas() {
NAMICS_DBG("LGrad1 computeLambda's " << endl);
	Real r, VL, LS;
	Real rlow, rhigh;

	if (fcc_sites){
		if (geometry=="cylindrical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=PIE*(pow(r,2)-pow(r-1,2));
				lambda1[i]=2.0*PIE*r/L[i]/3.0;
				lambda_1[i]=2.0*PIE*(r-1)/L[i]/3.0;
				lambda0[i]=1.0/3.0;
			}
		}
		if (geometry=="spherical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=4.0/3.0*PIE*(pow(r,3)-pow(r-1,3));
				fcc_lambda1[i]=4.0*PIE*pow(r,2)/L[i]/3.0;
				fcc_lambda_1[i]=4.0*PIE*pow(r-1,2)/L[i]/3.0;
				fcc_lambda0[i]=1.0-fcc_lambda1[i]-fcc_lambda_1[i];
			}
		}
	}

	if (fjc==1) {
		if (geometry=="cylindrical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=PIE*(pow(r,2)-pow(r-1,2));
				lambda1[i]=2.0*PIE*r/L[i]*lambda;
				lambda_1[i]=2.0*PIE*(r-1)/L[i]*lambda;
				lambda0[i]=1.0-2.0*lambda;
			}
		}
		if (geometry=="spherical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=4.0/3.0*PIE*(pow(r,3)-pow(r-1,3));
				lambda1[i]=4.0*PIE*pow(r,2)/L[i]*lambda;
				lambda_1[i]=4.0*PIE*pow(r-1,2)/L[i]*lambda;
				lambda0[i]=1.0-lambda1[i]-lambda_1[i];

			}
		}
		if (Markov==2) {
			for (int i=1; i<MX+1; i++) {
				l1[i]=lambda1[i]/lambda; l11[i]=1.0-l1[i];
				l_1[i]=lambda_1[i]/lambda; l_11[i]=1.0-l_1[i];
			}
		}

	}

	if (fjc>1) {
		if (geometry == "cylindrical") {
			for (int i = fjc; i < M - fjc; i++) {
				r = offset_first_layer+1.0*(i-fjc+1.0)/fjc;
				rlow = r - 0.5;
				rhigh = r + 0.5;
				L[i] = PIE * (2.0 * r) / fjc;
				VL = L[i] / PIE * fjc;
				if ((rlow - r) * 2 + r > 0.0)
					LAMBDA[i] += 1.0/(1.0*FJC-1.0)*rlow/VL;
				if ((rhigh - r) * 2 + r < 1.0*MX/fjc)
					LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
				else {
					if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2 * rhigh-r-1.0*MX/fjc < 0.001) {
						LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
					}
					for (int j = 1; j <= fjc; j++) {
						if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
							LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*(rhigh-1.0*j/fjc)/VL;
						}
					}
				}
				for (int j = 1; j < fjc; j++) {
					rlow += 0.5/(fjc);
					rhigh -= 0.5/(fjc);
					if ((rlow-r)*2+r > 0.0)
						LAMBDA[i+j*M] += 1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
					if ((rhigh-r)*2+r < offset_first_layer+1.0*MX/fjc)
						LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
					else {
						if (2 * rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
							LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
						}
						for (int k = 1; k <= fjc; k++) {
							if (2 * rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc<1.01*k/fjc) {
								LAMBDA[i + (FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
							}
						}
					}
				}
				LS = 0;
				for (int j = 0; j < FJC; j++)
					LS += LAMBDA[i+j*M];
				LAMBDA[i+(FJC/2)*M] += 1.0 - LS;
			}
		}

		if (geometry == "spherical") {
			for (int i = fjc; i < M - fjc; i++) {
				r = offset_first_layer+1.0*(1.0*i-1.0*fjc+1.0)/fjc;
				rlow = r-0.5;
				rhigh = r+0.5;
				L[i] = PIE*4.0/3.0*(rhigh*rhigh*rhigh-rlow*rlow*rlow)/fjc;
				VL = L[i] / PIE * fjc;
				if ((rlow-r)*2+r > 0.0)
					LAMBDA[i] += 0.5/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
				if ((rhigh -r)*2+r < 1.0*MX/fjc)
					LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
				else {
					if (2*rhigh-r-1.0*MX/fjc>-0.001 && 2*rhigh-r-1.0*MX/fjc<0.001) {
						LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
					}
					for (int j = 1; j <= fjc; j++) {
						if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
							LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*(rhigh-1.0*j/fjc)*(rhigh-1.0*j/fjc)/VL;
						}
					}
				}
				for (int j = 1; j < fjc; j++) {
					rlow += 0.5/(fjc);
					rhigh -= 0.5/(fjc);
					if ((rlow-r)*2+r > 0.0)
						LAMBDA[i+j*M] += 1.0/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
					if ((rhigh - r) * 2 + r < offset_first_layer + 1.0*MX/fjc)
						LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
					else {
						if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
							LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
						}
						for (int k = 1; k <= fjc; k++) {
							if (2*rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*k/fjc) {
								LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*(rhigh-1.0*k/fjc)*(rhigh-1.0*k/fjc)/VL;
							}
						}
					}
				}
				LS = 0;
				for (int j = 0; j < FJC; j++)
					LS += LAMBDA[i+j*M];
				LAMBDA[i+(FJC/2)*M] += 1.0-LS;
			}
		}
		if (Markov==2) { //planar gemaakt voor debugging....
			for (int i = fjc; i < M - fjc; i++) {
				LABDA[i]=LAMBDA[i]*8.0;
				LABDA_1[i]=1.0-LABDA[i];
				for (int j=1; j<FJC-1; j++) {
					LABDA[i+j*M]=LAMBDA[i+j*M]*4.0;
					LABDA_1[i+j*M]=1.0-LABDA[i+j*M];
				}
				LABDA[i+(FJC-1)*M]=LAMBDA[i+(FJC-1)*M]*8.0;
				LABDA_1[i+(FJC-1)*M]=1.0-LABDA[i+(FJC-1)*M];
			}
		}
	}
}

bool LGrad1::PutM() {
NAMICS_DBG("PutM in LGrad1 " << endl);	bool success=true;
	JX=1; JY=0; JZ=0; M=MX+2*fjc;
	if (geometry=="planar") {volume = MX/fjc; }
	if (geometry=="spherical") {volume = 4.0/3.0*PIE*(pow(MX+offset_first_layer,3)-pow(offset_first_layer,3))/fjc/fjc/fjc;}
	if (geometry=="cylindrical") {volume = PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2))/fjc/fjc;}

	Accesible_volume=volume;
	return success;
}

void LGrad1::TimesL(Real* X){
NAMICS_DBG("TimesL in LGrad1 " << endl); if (geometry!="planar") for (int __i = 0; __i < (M); ++__i) (X)[__i] = (X)[__i] * (L)[__i];
}

void LGrad1::DivL(Real* X){
NAMICS_DBG("DivL in LGrad1 " << endl); if (geometry!="planar") for (int __i = 0; __i < (M); ++__i) (X)[__i] = ((L)[__i] != 0) ? ((X)[__i] / (L)[__i]) : 0;
}

Real LGrad1:: Moment(Real* X,Real Xb, int n) {
NAMICS_DBG("Moment in LGrad1 " << endl);	Real Result=0;
	Real cor;
	remove_bounds(X);
	for (int i = fjc; i<M; i++) {
		cor = (i-fjc+0.5)/fjc;
		Result += pow(cor,n)*(X[i]-Xb)*L[i];
	}
	return Result/fjc;
}

Real LGrad1::MomentPlanar(Real* X,int n,Real Z0){
	(void)Z0;
	(void)n;
	(void)X;
	cout <<"MomentPlanar not implemented; kJ0 or kbar may be wrong. " << endl;
	return 0;
}

Real LGrad1::WeightedSum(Real* X){
NAMICS_DBG("weighted sum in LGrad1 " << endl);	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (X)[__i]; sum/=fjc;
	} else {
		(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (X)[__i] * (L)[__i];
	}
	return sum;
}

void LGrad1::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG(" Side in LGrad1 " << endl);
	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0);//set_bounds(X);
	int j, kk;

	if (fcc_sites) {
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i] * (fcc_lambda0)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i] * (fcc_lambda_1+1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i] * (fcc_lambda1)[__i];

	} else {
		if (fjc==1) {
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i] * (lambda0)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i] * (lambda_1+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i] * (lambda1)[__i];
		} else {

			for (j = 0; j < FJC/2; j++) {
				kk = (FJC-1)/2-j;
				for (int __i = 0; __i < (M-kk); ++__i) (X_side+kk)[__i] += (X)[__i] * (LAMBDA+j*M+kk)[__i];
				for (int __i = 0; __i < (M-kk); ++__i) (X_side)[__i] += (X+kk)[__i] * (LAMBDA+(FJC-j-1)*M)[__i];
			}
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i] * (LAMBDA+(FJC-1)/2*M)[__i];

		}
	}
}



void LGrad1::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	Real *gs=G+FJC*M*(s_to);
	Real *gs_1=G+FJC*M*(s_from);
	Real *g =G1;

	std::fill_n(gs, M*FJC, 0);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	if (fjc==1) {
		Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M;
		Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M;

		if (lattice_type==hexagonal) {

			int a,b; Real c;
			for (int p=0; p<FJC; p++){
				a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
				for (int q=0; q<FJC; q++) {
					c=P[abs(-p+q)];
					if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
					if (a>0) {
						for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] = (l_1+a)[__i] * (gs_1+q*M+b)[__i];
				  		for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] += (l_11+a)[__i] * (gs_1+(FJC-1-q)*M+a)[__i];
					}
					if (b>0) {
						for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] = (l1+a)[__i] * (gs_1+q*M+b)[__i];
				  		for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] += (l11+a)[__i] * (gs_1+(FJC-1-q)*M+a)[__i];
					}
					if (a+b>0){
				  		if (c!=0) for (int __i = 0; __i < (M-a-b); ++__i) (gs+p*M+a)[__i] += (c) * (H+b)[__i];
					} else {
						if (c!=0) for (int __i = 0; __i < (M); ++__i) (gs+p*M)[__i] += (c) * (gs_1+q*M)[__i];
					}
				}
			}
			for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		} else {
			for (int __i = 0; __i < (M-1); ++__i) (H)[__i] = (l_1 +1)[__i] * (gz0)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (H)[__i] += (l_11+1)[__i] * (gz2+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx0+1)[__i] += (P[0]) * (H)[__i];

				//LReflect(H,gz1,gz1);
			for (int __i = 0; __i < (M-1); ++__i) (H)[__i] = (l_1 +1)[__i] * (gz1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (H)[__i] += (l_11+1)[__i] * (gz1+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx0+1)[__i] += (4*P[1]) * (H)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (2*P[1]+P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (P[1]) * (gz2)[__i];

				//UReflect(H,gz1,gz1);
			for (int __i = 0; __i < (M-1); ++__i) (H+1)[__i] = (l1)[__i] * (gz1+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (H+1)[__i] += (l11)[__i] * (gz1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx2)[__i] += (4*P[1]) * (H+1)[__i];
				//UReflect(H,gz2,gz0);
			for (int __i = 0; __i < (M-1); ++__i) (H+1)[__i] = (l1)[__i] * (gz2+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (H+1)[__i] += (l11)[__i] * (gz0)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx2)[__i] += (P[0]) * (H+1)[__i];
			for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		int a,b; Real c;
		for (int p=0; p<FJC; p++){
			a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
			for (int q=0; q<FJC; q++) {
				c=P[abs(-p+q)];
				if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
				for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] = (LABDA+p*M+a)[__i] * (gs_1+q*M+b)[__i];
				for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] += (LABDA_1+p*M+a)[__i] * (gs_1+(FJC-1-q)*M+a)[__i];
				if (c!=0) for (int __i = 0; __i < (M-a-b); ++__i) (gs+p*M+a)[__i] += (c) * (H+b)[__i];
			}
		}
		for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
	}
}

void LGrad1::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	Real *gs=G+FJC*M*(s_to);
	Real *gs_1=G+FJC*M*(s_from);
	Real *g =G1;

	std::fill_n(gs, M*FJC, 0);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	if (fjc==1) {
		Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M;
		Real *gx0=gs,   *gx1=gs+M,   *gx2=gs+2*M;
		if (lattice_type==hexagonal) {
			int a,b; Real c;

			for (int q=FJC-1; q>-1; q--){
				a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
				if (a>0) {
					for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] = (l_1+a)[__i] * (gs_1+q*M+b)[__i];
					for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] += (l_11+a)[__i] * (gs_1+(FJC-1-q)*M+a)[__i];
				}
				if (b>0) {
					for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] = (l1+a)[__i] * (gs_1+q*M+b)[__i];
					for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] += (l11+a)[__i] * (gs_1+(FJC-1-q)*M+a)[__i];
				}
				for (int p=FJC-1; p>-1; p--) {
					c=P[abs(-p+q)];
					if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
					if (a+b>0) {
						if (c!=0) for (int __i = 0; __i < (M-a-b); ++__i) (gs+p*M+a)[__i] += (c) * (H+b)[__i];
					} else {
						if (c!=0) for (int __i = 0; __i < (M); ++__i) (gs+p*M)[__i] += (c) * (gs_1+q*M)[__i];
					}
				}
			}
			for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		} else {
			for (int __i = 0; __i < (M-1); ++__i) (H)[__i] = (l_1 +1)[__i] * (gz2)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (H)[__i] += (l_11+1)[__i] * (gz0+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx1+1)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx2+1)[__i] += (P[0]) * (H)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx0)[__i] += (4*P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (2*P[1]+P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (4*P[1]) * (gz1)[__i];

			for (int __i = 0; __i < (M-1); ++__i) (H+1)[__i] = (l1)[__i] * (gz0+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (H+1)[__i] += (l11)[__i] * (gz2)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx0)[__i] += (P[0]) * (H+1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gx1)[__i] += (P[1]) * (H+1)[__i];
			for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		int a,b; Real c;

		for (int q=FJC-1; q>-1; q--){
			a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
			for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] = (LABDA+(FJC-1-q)*M+a)[__i] * (gs_1+q*M+b)[__i];
			for (int __i = 0; __i < (M-a-b); ++__i) (H+b)[__i] += (LABDA_1+(FJC-1-q)*M+a)[__i] * (gs_1+(FJC-1-q)*M+a)[__i];
			for (int p=FJC-1; p>-1; p--) {
				c=P[abs(-p+q)];
				if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
				if (c!=0) for (int __i = 0; __i < (M-a-b); ++__i) (gs+p*M+a)[__i] += (c) * (H+b)[__i];
			}
		}
		for (int k=0; k<FJC; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
	}
}


void LGrad1::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
NAMICS_DBG(" propagate in LGrad1 " << endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int kk;
	int j;
	std::fill_n(gs, M, 0); set_bounds(gs_1);

	if (fjc==1) {
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (gs_1)[__i] * (lambda0)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (gs_1)[__i] * (lambda_1+1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (gs_1+1)[__i] * (lambda1)[__i];
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


bool LGrad1::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
	(void)var_pos;
	(void)n_pos;
	(void)H_p;
NAMICS_DBG("ReadRange in LGrad1 " << endl);	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	In->split(range,';',set);

	coor.clear();
	block=true;
	In->split(set[0],',',coor);
	if (coor.size()!=1) {cout << "In mon " + seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes must come as a single coordinate 'x'" << endl; r[0]=0; success=false;}
	else r[0]=ParseInt(coor[0],-1) ;

	coor.clear(); In->split(set[1],',',coor);

	if (coor.size()!=1) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "' the coordinates must come as a single coordinate 'x'" << endl; r[3]=0; success=false;}
	else r[3]=ParseInt(coor[0],-1);
	if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}

	return success;
}

bool LGrad1::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
NAMICS_DBG("ReadRangeFile in LGrad1 " << endl);	if (fjc>1) {
		cout << "Rangefile is not implemented for FJC-choices >3; contact FL. " << endl;
		return false;
	}
	bool success=true;
	string content;
	vector<string> lines;
	vector<string> xyz;
	string input_stem = In->name;
	const size_t slash_pos = input_stem.find_last_of("/\\");
	const size_t dot_pos = input_stem.find_last_of('.');
	if (dot_pos != string::npos && (slash_pos == string::npos || dot_pos > slash_pos)) {
		input_stem.erase(dot_pos);
	}
	const string resolved_range_file = input_stem + "." + filename;

	int length;
	int length_xyz;
	int px,p_i,x;
	int i=0;
	if (!io::ReadSanitizedFile(resolved_range_file,content)) {
		success=false;
		return success;
	}

	In->split(content,'#',lines);
	length = lines.size();
	if (length == MX) { //expect to read 'mask file';
		if (n_pos==0) {
			for (i = 0 ; i < length ; ++i) {
				if (ParseInt(lines[i],0)==1) n_pos++;
			}
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any elements." << endl;}
		} else {
			p_i=0;
			for (x=1; x<MX+1; x++) {
				if (ParseInt(lines[x-1],0)==1) {H_p[p_i]=x; p_i++;}
			}
		}
	} else { //expect to read x only
		px=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=1) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'single coordinate' 'x' was not found. " << endl;  success = false;
				} else {
					px=ParseInt(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
				}
				H_p[i]=px;
				i++;
			}
		}
	}
	return success;
}

bool LGrad1::FillMask(Real* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
	(void)pz;
	(void)py;
	bool success=true;
	bool readfile=false;
	int length=0;
	int length_px = px.size();

	vector<string> lines;
	int p;
	if (px.size()==0) {
		readfile=true;
		string content;
		success=io::ReadSanitizedFile(filename,content);
		if (success) {
			In->split(content,'#',lines);
			length = lines.size();
		}
	}
	if (readfile) {
		if (MX!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x-direction" << endl;
		} else {
			for (int x=1; x<MX+1; x++) Mask[x]=ParseInt(lines[x],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			else Mask[fjc-1+px[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad1::CreateMASK(Real* H_MASK, int* r, int* H_P, int n_pos, bool block) {
NAMICS_DBG("CreateMask for LGrad1 " + name << endl);	bool success=true;
	std::fill_n(H_MASK, M, static_cast<Real>(0));
	// Build mask from either a block in r=[x1,y1,z1,x2,y2,z2] or list of indices in H_P.
	if (block) {
		// mark all x in [x1, x2].
		for (int x=r[0]; x<r[3]+1; x++) {
			H_MASK[x]=1;
		}
	} else {
		// mark n_pos linear indices from H_P.
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad1::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	if (geometry !="planar") {
		result = 0;
		for (int __i = 0; __i < M; ++__i) result += phi[__i] * L[__i];
	} else {
		if (fjc==1) {
			result = 0;
			for (int __i = 0; __i < M; ++__i) result += phi[__i];
		} else {
			result = 0;
			for (int __i = 0; __i < M; ++__i) result += phi[__i] * L[__i];
		}
	}
	return result/fjc;
}

void LGrad1::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	std::fill_n(EE, M, 0);
	Real Exmin,Explus;
	int x;
	int r;

	if (geometry=="cylindrical" ) {
		r=offset_first_layer*fjc;
		pf=pf*PIE;
		for (x=fjc; x<MX+fjc; x++) {
			r++;
			Exmin=psi[x]-psi[x-1];
			Exmin*=(r-1)*Exmin;
			Explus=psi[x]-psi[x+1];
			Explus*=(r)*Explus;
			EE[x]=pf*(Exmin+Explus)/L[x];
		}




	}
	if (geometry=="spherical" ) {
		pf=pf*PIE*2/fjc;
		r=offset_first_layer*fjc +1.0;
		Explus=r*(psi[fjc]-psi[fjc+1]);
		Explus *=Explus;
		EE[fjc]=pf*Explus/(L[fjc]);
		for (x=fjc+1; x<MX+fjc; x++) {
			r +=1.0;
			Exmin=Explus;
			Explus=r*(psi[x]-psi[x+1]);
			Explus *=Explus;
			EE[x]=pf*(Exmin+Explus)/L[x];
		}

	}
}


void LGrad1::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, Real* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x;
	Real a,b,c,a_,b_,c_;
	Real r;
	Real epsXplus, epsXmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	if (geometry=="cylindrical") {
		C=C/PIE;
		r=offset_first_layer*fjc;
		epsXplus=r*(eps[fjc-1]+eps[fjc]);
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			r++;
			epsXmin=epsXplus;
			epsXplus=r*(eps[x]+eps[x+1]);
			//X[x]=(epsXmin*a + C*q[x]*L[x] + epsXplus*c)/(epsXmin+epsXplus);
			if (x==fjc) a=psi[fjc-1]; else a=X[x-1]; //upwind
			X[x]=(epsXmin*a  +C*q[x]*L[x] + epsXplus*psi[x+1])/(epsXmin+epsXplus);
		 }
	}
	if (geometry=="spherical") {
		C=C/(2.0*PIE)*fjc;
		r=offset_first_layer*fjc;
		epsXplus=r*r*(eps[fjc-1]+eps[fjc]);
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			epsXmin=epsXplus;
			r++;
			epsXplus=r*r*(eps[x]+eps[x+1]);
			a=b; b=c; c=psi[x+1];
			X[x]=(epsXmin*a + C*q[x]*L[x] + epsXplus*c)/(epsXmin+epsXplus);
		 }
	}
	for (int __i = 0; __i < (M); ++__i) (g)[__i] = (g)[__i] - (X)[__i];
   } else { //fixedPsi0 is true
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (x=fjc; x<MX+fjc; x++) {
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 0) psi[x]=0.5*(a+c)+q[x]*C/eps[x];
	}

	if (geometry=="cylindrical") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 0) psi[x]+=(c-a)/(2.0*(offset_first_layer*fjc+x-fjc+0.5))*fjc;
		}
	}
	if (geometry=="spherial") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 0) psi[x]+=(c-a)/(offset_first_layer*fjc+x-fjc+0.5)*fjc;
		}
	}
	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc];a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 0) {
				psi[x]+=0.25*(c_-a_)*(c-a)/eps[x]*fjc*fjc;
			}
		}
	}
	for (x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 0) {
		g[x]-=psi[x];
	}
   }
}


void LGrad1::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, Real* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x;
	Real a,b,c,a_,b_,c_;

	Real C = -e*e/(eps0*k_BT*bond_length);
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (x=fjc; x<MX+fjc; x++) { //for all geometries
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 1) q[x] = -0.5*(a-2*b+c)*fjc*fjc*eps[x]/C;
	}

	if (geometry=="cylindrical") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 1) q[x]-=(c-a)/(2.0*(offset_first_layer*fjc+x-fjc+0.5))*fjc*eps[x]/C;
		}
	}
	if (geometry=="spherial") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 1) q[x]-=(c-a)/(offset_first_layer*fjc+x-fjc+0.5)*fjc*eps[x]/C;
		}
	}
	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc]; a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 1) q[x]-=0.25*(c_-a_)*(c-a)*fjc*fjc/C;
		}
	}
	for (x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 1) {
		g[x]=-q[x];
	}

}

void LGrad1::remove_bounds(Real *X){
NAMICS_DBG("remove_bounds in LGrad1 " << endl);	int k;
	if (fjc==1) {
		X[0]=0;
		X[MX+1]=0;
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=0;
			X[MX+fjc+k]=0;
		}
	}
}

void LGrad1::set_bounds(Real* X, Real* Y){
NAMICS_DBG("set_bounds in LGrad1 " << endl);	int k;
	if (fjc==1) {
		X[0]=Y[BX1];
		X[MX+1]=Y[BXM];
		Y[0]=X[BX1];
		Y[MX+1]=X[BXM];

	} else {
		for (k=0; k<fjc; k++) {
			X[k]=Y[B_X1[k]];
			X[MX+fjc+k]=Y[B_XM[k]];
			Y[k]=X[B_X1[k]];
			Y[MX+fjc+k]=X[B_XM[k]];
		}
	}
}

void LGrad1::set_bounds(Real* X){
NAMICS_DBG("set_bounds in LGrad1 " << endl);	int k=0;
	if (fjc==1) {
		X[0]=X[BX1];
		X[MX+1]=X[BXM];
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=X[B_X1[k]];
			X[MX+fjc+k]=X[B_XM[k]];
		}
	}
}

void LGrad1::set_M_bounds(Real* X){
NAMICS_DBG("set_M_bounds in LGrad1 " << endl); //set mirror bounds
	int k=0;
	if (fjc==1) {
		X[0]=X[1];
		X[MX+1]=X[MX];
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=X[2*fjc-k-1];
			X[MX+fjc+k]=X[MX+fjc-k-1];
		}
	}
}


void LGrad1::remove_bounds(int *X){
NAMICS_DBG("remove_bounds in LGrad1 " << endl);int k;
	if (fjc==1) {
		X[0]=0;
		X[MX+1]=0;
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=0;
			X[MX+fjc+k]=0;
		}
	}
}

void LGrad1::set_bounds(int* X){
NAMICS_DBG("set_bounds in LGrad1 " << endl);	int k=0;
	if (fjc==1) {
		X[0]=X[BX1];
		X[MX+1]=X[BXM];
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=X[B_X1[k]];
			X[MX+fjc+k]=X[B_XM[k]];
		}
	}
}

Real LGrad1::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	if (Markov==2) {
		GN=WeightedSum(G);
		for (int k=1; k<FJC-1; k++) {
			if (lattice_type == hexagonal) GN += 2.0*WeightedSum(G+k*M); else GN +=4.0*WeightedSum(G+k*M);
		}
		GN+=WeightedSum(G+(FJC-1)*M);
		if (lattice_type == hexagonal) GN /= 4.0*fjc; else GN /= 6.0;
	} else GN=WeightedSum(G);
	return GN;
}

void LGrad1::AddPhiS(Real* phi,Real* Gf,Real* Gb,int Markov, int M){
	NAMICS_DBG("AddPhiS_markov " << endl);	if (Markov==2) {


		if (lattice_type ==hexagonal) {
			if (fjc==1) {
				Real C1=1.0/4.0;
				Real C2=2.0/4.0;
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf)[__i] * (Gb)[__i];
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C2) * (Gf+1*M)[__i] * (Gb+1*M)[__i];
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf+2*M)[__i] * (Gb+2*M)[__i];
			} else {
				Real C1=0.5/(FJC-1.0);
				Real C2=1.0/(FJC-1.0);
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf)[__i] * (Gb)[__i];
				for (int k=1; k<FJC-1; k++) for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C2) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf+(FJC-1)*M)[__i] * (Gb+(FJC-1)*M)[__i];
			}
		} else { //markov=2 cubic fjc=1
			Real C1=1.0/6.0;
			Real C2=4.0/6.0;
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf)[__i] * (Gb)[__i];
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C2) * (Gf+1*M)[__i] * (Gb+1*M)[__i];
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf+2*M)[__i] * (Gb+2*M)[__i];
		}
	} else {
		for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (Gf)[__i] * (Gb)[__i];
	}
}

void LGrad1::AddPhiS(Real* phi,Real* Gf,Real* Gb, Real degeneracy, int Markov, int M){
NAMICS_DBG("AddPhiS_degeneracy markov " << endl);	if (Markov==2) {
		if (lattice_type ==hexagonal) {
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy*0.5/(FJC-1.0)) * (Gf)[__i] * (Gb)[__i];
			for (int k=1; k<FJC-1; k++) for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/(FJC-1.0)) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy*0.5/(FJC-1.0)) * (Gf+(FJC-1)*M)[__i] * (Gb+(FJC-1)*M)[__i];
		} else {
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/6.0) * (Gf)[__i] * (Gb)[__i];
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy*4.0/6.0) * (Gf+1*M)[__i] * (Gb+1*M)[__i];
			for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/6.0) * (Gf+2*M)[__i] * (Gb+2*M)[__i];
		}
	} else {
		for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy) * (Gf)[__i] * (Gb)[__i];
	}
}

void LGrad1::Initiate(Real* G,Real* Gz,int Markov, int M){
	if (Markov==2) {
		for (int k=0; k<FJC; k++) std::copy_n(Gz, M, G+k*M);
	} else {
		std::copy_n(Gz, M, G);
	}
}

void LGrad1::Terminate(Real* Gz ,Real* G, int Markov, int M){
NAMICS_DBG("LGrad1::Terminate " << endl);	Real one=1.0;
	if (Markov==2) {
		std::fill_n(Gz, M, 0);
		if (lattice_type == simple_cubic) {
			for (int __i = 0; __i < M; ++__i) Gz[__i] += (G + M)[__i];
			for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (4.0*one);
			for (int __i = 0; __i < M; ++__i) Gz[__i] += G[__i];
			for (int __i = 0; __i < M; ++__i) Gz[__i] += (G + 2 * M)[__i];
			for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (1.0/6.0*one);
		} else {
			for (int __i = 0; __i < M; ++__i) Gz[__i] += (G + M)[__i];
			for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (2.0*one);
			for (int __i = 0; __i < M; ++__i) Gz[__i] += G[__i];
			for (int __i = 0; __i < M; ++__i) Gz[__i] += (G + 2 * M)[__i];
			for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (1.0/4.0*one);

		}
	} else {
		std::copy_n(G, M, Gz);
	}
}

bool LGrad1:: PutMask(Real* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	(void)R;
	(void)pz;
	(void)py;
	(void)px;
	(void)MASK;
	bool success=false;
	cout <<"PutMask does not make sence in 1 gradient system " << endl;
	return success;
}
