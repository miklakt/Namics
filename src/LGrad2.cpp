#include <iostream>
#include <string>
#include "lattice.h"
#include "LGrad2.h"

LGrad2::LGrad2(const Input& In_,const string& name_): Lattice(In_,name_) {}

LGrad2::~LGrad2() {
NAMICS_DBG("LGrad2 destructor " << endl);}


void LGrad2:: ComputeLambdas() {
	Real r, VL, LS;
	Real rlow, rhigh;


	//		L[P(x,y)]=PIE*(pow(r,2)-pow(r-1,2));
	//		fcc_lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]/3.0;
	//		fcc_lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]/3.0;
	//		fcc_lambda0[P(x,y)]=1.0-2.0/3.0;
	//	}
	//
	//}
	if (fjc==1) {
		for (int x=1; x<MX+1; x++)
		for (int y=1; y<MY+1; y++) {
			r=offset_first_layer + 1.0*x;
			L[P(x,y)]=PIE*(pow(r,2)-pow(r-1,2));
			lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]*lambda;
			lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]*lambda;
			lambda0[P(x,y)]=1.0-2.0*lambda;
			if (fcc_sites) {
				fcc_lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]/3.0;
				fcc_lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]/3.0;
				fcc_lambda0[P(x,y)]=1.0-2.0/3.0;
			}
		}
		if (Markov ==2) {
			for (int i=0; i<M; i++) {
				l1[i]=lambda1[i]/lambda; l11[i]=1.0-l1[i];
				l_1[i]=lambda_1[i]/lambda; l_11[i]=1.0-l_1[i];
			}
		}

	}
	if (fjc>1) {
		for (int y=fjc; y<MY+fjc; y++) {
			for (int x = fjc; x < MX+fjc; x++) {
				r = offset_first_layer+1.0*(x-fjc+1.0)/fjc;
				rlow = r - 0.5;
				rhigh = r + 0.5;
				L[P(x,y)] = PIE * (2.0 * r) / fjc;
				VL = L[P(x,y)] / PIE * fjc;
				if ((rlow - r) * 2 + r > 0.0) {
					LAMBDA[P(x,y)] += 1.0/(1.0*FJC-1.0)*rlow/VL;
				}
				if ((rhigh - r) * 2 + r < 1.0*MX/fjc) {
					LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
				} else {
					if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2 * rhigh-r-1.0*MX/fjc < 0.001) {
						LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
					}
					for (int j = 1; j <= fjc; j++) {
						if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
							LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*(rhigh-1.0*j/fjc)/VL;
						}
					}								}
				for (int j = 1; j < fjc; j++) {
					rlow += 0.5/(fjc);
					rhigh -= 0.5/(fjc);
					if ((rlow-r)*2+r > 0.0)
					LAMBDA[P(x,y)+j*M] += 1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
					if ((rhigh-r)*2+r < offset_first_layer+1.0*MX/fjc)
					LAMBDA[P(x,y)+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
					else {
						if (2 * rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
							LAMBDA[P(x,y)+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
						}
						for (int k = 1; k <= fjc; k++) {
							if (2 * rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc<1.01*k/fjc) {
								LAMBDA[P(x,y) + (FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
							}
						}
					}
				}
				LS = 0;
				for (int j = 0; j < FJC; j++)
				LS += LAMBDA[P(x,y)+j*M];
				LAMBDA[P(x,y)+(FJC/2)*M] += 1.0 - LS;
			}
		}
	}
		//LAMBDA[P(x,y)+k*M]=LAMBDA[P(2*MX-x+1,y)+(FJC-k-1)*M];
}

bool LGrad2::PutM() {
NAMICS_DBG("PutM in LGrad2 " << endl);	bool success=true;

	if (geometry=="cylindrical")
		volume = MY*PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));
	else volume = MX*MY;
	JX=MY+2*fjc; JY=1; JZ=0; M=(MX+2*fjc)*(MY+2*fjc);

	Accesible_volume=volume;
	return success;
}

void LGrad2::TimesL(Real* X){
NAMICS_DBG("TimesL in LGrad2 " << endl); if (geometry!="planar") for (int __i = 0; __i < (M); ++__i) (X)[__i] = (X)[__i] * (L)[__i];
}

void LGrad2::DivL(Real* X){
NAMICS_DBG("DivL in LGrad2 " << endl); if (geometry!="planar") for (int __i = 0; __i < (M); ++__i) (X)[__i] = ((L)[__i] != 0) ? ((X)[__i] / (L)[__i]) : 0;
}

Real LGrad2:: Moment(Real* X,Real Xb, int n) {
NAMICS_DBG("Moment in LGrad2 " << endl);	Real Result=0;
	int x,y;
	int zz;
	Real Nz;
	if (fjc==1) {
		zz=0;
		for (y=1; y<=MY ; y++) {
			Nz=0;
			for (x=1; x<=MX; x++) {
				if (X[P(x,y)]>0) Nz+=(X[P(x,y)]-Xb)*L[P(x,y)];
			}
			if (Nz>0) zz++;
			if (zz>0) Result+= pow(zz,n)*Nz;
		}
	} else {
	}
	return Result/fjc;
}

Real LGrad2::WeightedSum(Real* X){
NAMICS_DBG("weighted sum in LGrad2 " << endl);	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (X)[__i]; sum = sum/(fjc*fjc);
	} else	{
		(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (X)[__i] * (L)[__i];
	}
	return sum;
}

void LGrad2::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG(" Side in LGrad2 " << endl);	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0);//set_bounds(X);

	if (fcc_sites) {
		Real C1=1.0/3.0;
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (C1) * (X)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i] * (fcc_lambda_1+JX)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i] * (fcc_lambda1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (C1) * (X)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (C1) * (X+1)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (X)[__i] * (fcc_lambda_1+JX+1)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX)[__i] += (X+1)[__i] * (fcc_lambda_1+JX)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+1)[__i] += (X+JX)[__i] * (fcc_lambda1+1)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JX+1)[__i] * (fcc_lambda1)[__i];
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C1);

	} else {
		if (fjc==1) {
			if (lattice_type ==simple_cubic) {
				Real C1=4.0/6.0;
				Real C2=1.0/6.0;
				Real C3=4.0;
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (C1) * (X)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i] * (lambda_1+JX)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i] * (lambda1)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (C2) * (X)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (C2) * (X+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C3);
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (X)[__i] * (lambda_1+JX+1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX)[__i] += (X+1)[__i] * (lambda_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+1)[__i] += (X+JX)[__i] * (lambda1+1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JX+1)[__i] * (lambda1)[__i];
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C2);
			} else {
				Real C1=2.0/4.0;
				Real C2=1.0/4.0;
				Real C3=2.0;
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (C1) * (X)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i] * (lambda_1+JX)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i] * (lambda1)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (C2) * (X)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (C2) * (X+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C3);
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (X)[__i] * (lambda_1+JX+1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX)[__i] += (X+1)[__i] * (lambda_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+1)[__i] += (X+JX)[__i] * (lambda1+1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JX+1)[__i] * (lambda1)[__i];
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C2);
			}
		} else { //fjc>1
			for (int block=0; block<2; block++) {
				int a,b;
				int bk;
				Real C=1.0/(2*(fjc+2));
				for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) {
					bk=0; a=0; b=0;
					if (y==-fjc || y==fjc) bk++;
					if (bk==block) {
						if (x<0) a =-x*JX;  else b=x*JX;
						if (y<0) a -=y*JY;  else b+=y*JY;
						for (int __i = 0; __i < (M-a-b); ++__i) (X_side+a)[__i] += (X+b)[__i] * (LAMBDA+(fjc+x)*M+a)[__i];
					}
				}
				if (block !=1) for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (2.0); else for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);
			}
		}
	}
}


void LGrad2::LReflect(Real *H, Real *P, Real *Q) {
	for (int __i = 0; __i < (M-JX); ++__i) (H)[__i] = (l_1+JX)[__i] * (P)[__i];
	for (int __i = 0; __i < (M-JX); ++__i) (H)[__i] += (l_11+JX)[__i] * (Q+JX)[__i];
}

void LGrad2::UReflect(Real *H, Real *P, Real *Q) {
	for (int __i = 0; __i < (M-JX); ++__i) (H+JX)[__i] = (l1)[__i] * (P+JX)[__i];
	for (int __i = 0; __i < (M-JX); ++__i) (H+JX)[__i] += (l11)[__i] * (Q)[__i];
}


void LGrad2::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {
		if (lattice_type == hexagonal) {
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs,   *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0);set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			LReflect(H,gz0,gz11);for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (H)[__i]; //0 and 1 are equivalent
			LReflect(H,gz1,gz10);for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (H)[__i]; //10 and 11 are equivalent
			LReflect(H,gz2,gz9);for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (H)[__i]; //3 and 4 equivalent
								for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz3)[__i]; //7 and 8

								for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz4)[__i];
			//LReflect(H,gz5,gz6);
								for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz5)[__i];
			//LReflect(H,gz6,gz5);
								for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz6)[__i];
			//LReflect(H,gz7,gz4);
								for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz7)[__i];
			//LReflect(H,gz8,gz3);
								for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz8)[__i];


			//UReflect(H,gz3,gz8);
								for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz3+JX)[__i];
			//UReflect(H,gz4,gz7);
								for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz4+JX)[__i];
			//UReflect(H,gz5,gz6);
								for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz5+JX)[__i];
			//UReflect(H,gz6,gz5);
								for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz6+JX)[__i];
			//UReflect(H,gz7,gz4);
								for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz7+JX)[__i];
			//UReflect(H,gz8,gz3);
								for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz8+JX)[__i];
			UReflect(H,gz9,gz2);for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[0]) * (H+JX)[__i];
			UReflect(H,gz10,gz1);for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[0]) * (H+JX)[__i];
			UReflect(H,gz11,gz0);for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[0]) * (H+JX)[__i];


			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz2,gz9,0);  set_bounds_y(gz3,gz8,0); set_bounds_y(gz4,gz7,0); set_bounds_y(gz0,gz11,0); set_bounds_y(gz1,gz10,0); set_bounds_y(gz5,gz6,0);

			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[0]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz11)[__i];


			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[0]) * (gz6+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[0]) * (gz7+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[0]) * (gz8+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz10+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz11+JY)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[0]) * (gz4)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz11)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[0]) * (gz7)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[0]) * (gz8)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz11)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			LReflect(H,gz0,gz11);for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[0]) * (H)[__i];
			LReflect(H,gz1,gz10);for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[0]) * (H)[__i];
			LReflect(H,gz2,gz9);for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[0]) * (H)[__i];
			//LReflect(H,gz3,gz8);
								for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz3)[__i];
			//LReflect(H,gz4,gz7);
								for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz4)[__i];
			//LReflect(H,gz5,gz6);
								for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz5)[__i];
			//LReflect(H,gz6,gz5);
								for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz6)[__i];
			//LReflect(H,gz7,gz4);
								for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz7)[__i];
			//LReflect(H,gz8,gz3);
								for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz8)[__i];

			//UReflect(H,gz3,gz8);
								for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz3+JX)[__i];
			//UReflect(H,gz4,gz7);
								for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz4+JX)[__i];
			//UReflect(H,gz5,gz6);
								for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz5+JX)[__i];
			//UReflect(H,gz6,gz5);
								for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz6+JX)[__i];
			//UReflect(H,gz7,gz4);
								for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz7+JX)[__i];
			//UReflect(H,gz8,gz3);
								for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz8+JX)[__i];
			UReflect(H,gz9,gz2);for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[0]) * (H+JX)[__i];
			UReflect(H,gz10,gz1);for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[0]) * (H+JX)[__i];
			UReflect(H,gz11,gz0);for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[0]) * (H+JX)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz2,gz9,0);  set_bounds_y(gz3,gz8,0); set_bounds_y(gz4,gz7,0);
			set_bounds_y(gz0,gz11,0); set_bounds_y(gz1,gz10,0); set_bounds_y(gz5,gz6,0);

			//LReflect(H,gz0,gz11);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz0+JY)[__i];
			//LReflect(H,gz1,gz10);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz1+JY)[__i];
			//LReflect(H,gz2,gz9);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz2+JY)[__i];
			//LReflect(H,gz3,gz8);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz3+JY)[__i];
			//LReflect(H,gz4,gz7);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz4+JY)[__i];
			//LReflect(H,gz5,gz6);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz5+JY)[__i];
			//LReflect(H,gz6,gz5);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz6+JY)[__i];
			//LReflect(H,gz7,gz4);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz7+JY)[__i];
			//LReflect(H,gz8,gz3);
								for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz8+JY)[__i];

			//UReflect(H,gz3,gz8);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz3)[__i];
			//UReflect(H,gz4,gz7);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz4)[__i];
			//UReflect(H,gz5,gz6);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz5)[__i];
			//UReflect(H,gz6,gz5);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz6)[__i];
			//UReflect(H,gz7,gz4);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz7)[__i];
			//UReflect(H,gz8,gz3);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz8)[__i];
			//UReflect(H,gz9,gz2);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz9)[__i];
			//UReflect(H,gz10,gz1);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz10)[__i];
			//UReflect(H,gz11,gz0);
								for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz11)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[0]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz11)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[0]) * (gz6+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[0]) * (gz7+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[0]) * (gz8+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz10+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz11+JY)[__i];

			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];



		} else {//simple _cubic should work
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M,*gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M,*gx4=gs+4*M;
			Real *g=G1;

			std::fill_n(gs, 5*M, 0);
			LReflect(H,gz0,gz4); for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (2*P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz3)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (2*P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz4)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[1]) * (gz4)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[1]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (2*P[1]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[0]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[1]) * (gz4+JY)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (2*P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (gz3+JX)[__i];
			UReflect(H,gz4,gz0); for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[0]) * (H+JX)[__i];

			for (int k=0; k<5; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		if (lattice_type==hexagonal) {
			cout <<"stencil_full, hexagonal, Markov 2, cyl coordinates not implemented" << endl;
		} else {
			cout <<"stencil_full, simple_cubic, Markov 2, cyl coordinates not implemented" << endl;
		}
	}
}

void LGrad2::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {

	if (!stencil_full) {
		if (lattice_type==hexagonal) {
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			for (int k=0; k<12; k++) remove_bounds(gs_1+k*M);

			set_bounds_x(gz0,gz11,0);
			LReflect(H,gz11,gz0);
			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx9+JX)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10+JX)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11+JX)[__i] += (P[0]) * (H)[__i];

			UReflect(H,gz0,gz11);
			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8)[__i] += (P[1]) * (H+JX)[__i];

			remove_bounds(gz0);remove_bounds(gz11);
			set_bounds_y(gz3,gz8,0);

			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx6+JY)[__i] += (P[0]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7+JY)[__i] += (P[0]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8+JY)[__i] += (P[0]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11+JY)[__i] += (P[1]) * (gz8)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx0)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[0]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[0]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5)[__i] += (P[0]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11)[__i] += (P[1]) * (gz3+JY)[__i];

			remove_bounds(gz3);remove_bounds(gz8);

			for (int __i = 0; __i < (M); ++__i) (gx0)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx7)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx8)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx9)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx10)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx11)[__i] += (P[1]) * (gz6)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx0)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx3)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx4)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx9)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx10)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx11)[__i] += (P[1]) * (gz5)[__i];

			remove_bounds(gz5); remove_bounds(gz6);
			set_bounds_x(gz1,gz10,0);

			LReflect(H,gz10,gz1);

			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx9+JX)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10+JX)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11+JX)[__i] += (P[0]) * (H)[__i];

			UReflect(H,gz1,gz10);

			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8)[__i] += (P[1]) * (H+JX)[__i];

			remove_bounds(gz1);remove_bounds(gz10);
			//set_bounds_y(gz2,gz9,0);
			//LReflect(H,gz9,gz2);
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[1]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5)[__i] += (P[1]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx6)[__i] += (P[1]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9)[__i] += (P[0]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[0]) * (H+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11)[__i] += (P[0]) * (H+JY)[__i];

			//UReflect(H,gz2,gz9);
			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[0]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5+JY)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx6+JY)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7+JY)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8+JY)[__i] += (P[1]) * (H)[__i];

			remove_bounds(gz2);remove_bounds(gz9);

			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx6+JY)[__i] += (P[0]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7+JY)[__i] += (P[0]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8+JY)[__i] += (P[0]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10+JY)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11+JY)[__i] += (P[1]) * (gz7)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx0)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11)[__i] += (P[1]) * (gz4+JY)[__i];

			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];


		} else { //simple_cubic should work
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
			Real *g=G1;

			std::fill_n(gs, 5*M, 0);

			LReflect(H,gz4,gz0);
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (H)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[0]) * (H)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz3)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx0)[__i] += (2*P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx1)[__i] += (2*P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx2)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx3)[__i] += (2*P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx4)[__i] += (2*P[1]) * (gz2)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx0)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1)[__i] += (P[0]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (gz1+JY)[__i];

			UReflect(H,gz0,gz4);
			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[1]) * (H+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (H+JX)[__i];

			for (int k=0; k<5; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		if (lattice_type==hexagonal) {
			cout <<"stencil_full, cyl coordinates, hexagonal, Markov 2 not implemented" << endl;
		} else {
			cout <<"stencil_full, cyl coordinates, simple_cubic, Markov 2 not implemented" << endl;
		}
	}
}

void LGrad2::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
NAMICS_DBG(" propagate in LGrad2 " << endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	std::fill_n(gs, M, 0); set_bounds(gs_1);
	if (fjc==1) {
		if (lattice_type==simple_cubic) {
			Real C1=4.0/6.0;
			Real C2=1.0/6.0;
			Real C3=4.0;
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (C1) * (gs_1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i] * (lambda_1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i] * (lambda1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (C2) * (gs_1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (C2) * (gs_1+1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C3);
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX+1)[__i] += (gs_1)[__i] * (lambda_1+JX+1)[__i];
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX)[__i] += (gs_1+1)[__i] * (lambda_1+JX)[__i];
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs+1)[__i] += (gs_1+JX)[__i] * (lambda1+1)[__i];
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs)[__i] += (gs_1+JX+1)[__i] * (lambda1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C2);
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
		} else { //9 point stencil; hexagonal
			Real C1=0.5;
			Real C2=0.25;
			Real C3=2.0;
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (C1) * (gs_1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i] * (lambda_1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i] * (lambda1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (C2) * (gs_1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (C2) * (gs_1+1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C3);
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX+1)[__i] += (gs_1)[__i] * (lambda_1+JX+1)[__i];
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX)[__i] += (gs_1+1)[__i] * (lambda_1+JX)[__i];
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs+1)[__i] += (gs_1+JX)[__i] * (lambda1+1)[__i];
			for (int __i = 0; __i < (M-JX-1); ++__i) (gs)[__i] += (gs_1+JX+1)[__i] * (lambda1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C2);
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
		}
	} else { //fjc>1
		for (int block=0; block<2; block++) {
			int a,b;
			int bk;
			Real C= 1.0/(2*(fjc+2));
			for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) {
				bk=0; a=0; b=0;
				if (y==-fjc || y==fjc) bk++;
				if (bk==block) {
					if (x<0) a =-x*JX;  else b=x*JX;
					if (y<0) a -=y*JY;  else b+=y*JY;
					for (int __i = 0; __i < (M-a-b); ++__i) (gs+a)[__i] += (gs_1+b)[__i] * (LAMBDA+(fjc+x)*M+a)[__i];
				}
			}
			if (block !=1) for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (2.0); else for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C);
		}
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];


	}
}


bool LGrad2::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
	(void)var_pos;
	(void)n_pos;
	(void)H_p;
NAMICS_DBG("ReadRange in LGrad2 " << endl);	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	In->split(range,';',set);
	coor.clear();
	block=true; In->split(set[0],',',coor);
	if (coor.size()!=2) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of two: 'x,y'" << endl; success=false;}
	else {
		r[0]=ParseInt(coor[0],0);
		r[1]=ParseInt(coor[1],0);
	}
	coor.clear(); In->split(set[1],',',coor);

	if (coor.size()!=2) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of two: 'x,y'" << endl; success=false;}
	else {
		r[3]=ParseInt(coor[0],0);
		r[4]=ParseInt(coor[1],0);
	}
	if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}

	return success;
}

bool LGrad2::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
NAMICS_DBG("ReadRangeFile in LGrad2 " << endl);	if (fjc>1) {
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
	int px,py,p_i,x,y;
	int i=0;
	if (!io::ReadSanitizedFile(resolved_range_file,content)) {
		success=false;
		return success;
	}

	In->split(content,'#',lines);
	length = lines.size();


	if (length == MX*MY) { //expect to read 'mask file';
		i=0;
		if (n_pos==0) {
			while (i<length){
				if (ParseInt(lines[i],0)==1) n_pos++;
				i++;
			};
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
		} else {
			i=0; p_i=0;
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++)  {
				if (ParseInt(lines[i],0)==1) {H_p[p_i]=P(x,y); p_i++;}
				i++;
			}
		}
	} else { //expect to read x,y
		px=0,py=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=2) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'pair of coordinates' 'x,y' was not found. " << endl;  success = false;
				} else {
					px=ParseInt(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
					py=ParseInt(xyz[1],0);
					if (py < 1 || py > MY) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MY << endl; success =false;}
				}
				cout <<"reading px " << px << " and py " << py << endl;
				H_p[i]=P(px,py);
				i++;
			}
		}
	}
	return success;
}

bool LGrad2::FillMask(Real* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
	(void)pz;
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
		if (MX*MY!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x,y-directions" << endl;
		} else {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++) Mask[x*JX + y]=ParseInt(lines[x*JX+y],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			p=py[i]; if (p<1 || p>MY) {success=false; cout<<" y-value in delta_range out of bounds; " << endl; }
			if (success) Mask[P(px[i],py[i])]=1; //Mask[px[i]*JX + fjc-1+ py[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad2::CreateMASK(Real* H_MASK, int* r, int* H_P, int n_pos, bool block) {
NAMICS_DBG("CreateMask for LGrad2 " + name << endl);	bool success=true;
	std::fill_n(H_MASK, M, static_cast<Real>(0));
	// Build mask from either a block in r=[x1,y1,z1,x2,y2,z2] or list of indices in H_P.
	if (block) {
		// mark all (x,y) in [x1, x2] x [y1, y2].
		for (int x=r[0]; x<r[3]+1; x++)
		for (int y=r[1]; y<r[4]+1; y++)
			H_MASK[P(x,y)]=1;
	} else {
		// mark n_pos linear indices from H_P.
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad2::ComputeTheta(Real* phi) {
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
	return result;
}

void LGrad2::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	std::fill_n(EE, M, 0);
	Real Exmin,Explus,Eymin,Eyplus;
	int x,y,z;
	int r;

	pf = pf/2;
	r=offset_first_layer*fjc;
	for (x=fjc; x<MX+fjc; x++) {
		r++;
		for (y=fjc; y<MY+fjc; y++) {
			z=x*JX+y;
			Exmin=psi[z]-psi[z-JX];
			Exmin*=(r-0.5)*Exmin*2*PIE;
			Explus=psi[z]-psi[z+JX];
			Explus*=(r+0.5)*Explus*2*PIE;
			Eymin=(psi[z]-psi[z-1])*fjc;
			Eymin*=Eymin;
			Eyplus=(psi[z]-psi[z+1])*fjc;
			Eyplus*=Eyplus;
			EE[z]=pf*((Exmin+Explus)/L[z]+Eymin+Eyplus);
		}
	}
}


void LGrad2::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, Real* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y,i;

	Real r;
	Real epsXplus, epsXmin, epsYplus,epsYmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);
	Real a=0;
	Real b=0;
	if (!fixedPsi0) {
		C=C/2/fjc/fjc;
		r=offset_first_layer*fjc;
		for (x=fjc; x<MX+fjc; x++) {
			r++;
			for (y=fjc; y<MY+fjc; y++) {
				i=x*JX+y;
				epsXmin=2*PIE*(r-1)*(eps[i]+eps[i-JX])/L[i]*fjc*fjc;
				epsXplus=2*PIE*r*(eps[i]+eps[i+JX])/L[i]*fjc*fjc;
				epsYmin=eps[i]+eps[i-1];
				epsYplus=eps[i]+eps[i+1];
				if (x==fjc) a=psi[i-JX]; else a=X[i-JX];
				if (y==fjc) b=psi[i-1]; else b=X[i-1];
				X[i]= (C*q[i]+epsXmin*a+epsXplus*psi[i+JX]+epsYmin*b+epsYplus*psi[i+1])/
				(epsXmin+epsXplus+epsYmin+epsYplus);
				//X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/
				//(epsXmin+epsXplus+epsYmin+epsYplus);
			}
		}
		for (int __i = 0; __i < (M); ++__i) (g)[__i] = (g)[__i] - (X)[__i];
  	 } else { //fixedPsi0 is true
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]=0.25*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
			        +psi[x*JX+y-1]  +psi[x*JX+y+1])
				 +0.5*q[x*JX+y]*C/eps[x*JX+y];
			}
		}


		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++){
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]+=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc;
			}
		}
		if (grad_epsilon) {
			for (x=fjc; x<MX+fjc; x++) {
				for (y=fjc; y<MY+fjc; y++) {
					if (Mask[x*JX+y] == 0) {
						X[x*JX+y]+=0.25*(eps[(x+1)*JX+y]-eps[(x-1)*JX+y])*(psi[(x+1)*JX+y]-psi[(x-1)*JX+y]+
                                              eps[x*JX+y+1]  -eps[x*JX+y-1])  *(psi[x*JX+y+1]  -psi[x*JX+y-1])/
					           eps[x*JX+y]*fjc*fjc;
					}
				}
			}
		}
		for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)
		if (Mask[x*JX+y] == 0) {
			psi[x*JX+y]=X[x*JX+y];
			g[x*JX+y]-=psi[x*JX+y];
		}
   	}
}


void LGrad2::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, Real* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x,y;

	Real C = -e*e/(eps0*k_BT*bond_length);
	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++){ //for all geometries
			if (Mask[x*JX+y] == 1)
			q[x*JX+y]=-0.5*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
				        +psi[x*JX+y-1]  +psi[x*JX+y+1]
					 -4*psi[x*JX+y])*fjc*fjc*eps[x*JX+y]/C;
		}
	}

	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++){
			if (Mask[x*JX+y] == 1)
			q[x*JX+y]-=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc*eps[x]/C;
		}
	}
	if (grad_epsilon) {
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				if (Mask[x*JX+y] == 1)
				q[x*JX+y]-=0.25*(eps[(x+1)*JX+y]-eps[(x-1)*JX+y])*(psi[(x+1)*JX+y]-psi[(x-1)*JX+y]+
                                  		   eps[x*JX+y+1]  -eps[x*JX+y-1])  *(psi[x*JX+y+1]  -psi[x*JX+y-1])*fjc*fjc/C;
				}
		}
	}
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)
	if (Mask[x*JX+y] == 1) {
		g[x*JX+y]=-q[x*JX+y];
	}

}

void LGrad2::remove_bounds(Real *X){
NAMICS_DBG("remove_bounds in LGrad2 " << endl);	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = 0;
			X[x*JX+MY+1]=0;
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = 0;
			X[(MX+1)*JX+y]=0;
		}
		//corners
		for (x=0; x<1; x++) {
			X[x*JX+0] = 0;
			X[x*JX+MY+1]=0;
		}
		for (x=MX+1; x<MX+2; x++) {
			X[x*JX+0] = 0;
			X[x*JX+MY+1]=0;
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=0;
				X[x*JX+MY+fjc+k]=0;
			}
		}
		for (y=0; y<MY+2*fjc; y++) { //this will also remove the corners...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=0;
				X[(MX+fjc+k)*JX+y]=0;
			}
		}
	}
}

void LGrad2::set_bounds_x(Real* X, Real*Y, int shifty){
NAMICS_DBG("set_bounds_x XY in LGrad2 " << endl);
	if (BX1>BXM)  {
		//set_bounds_x(X,0); set_bounds_x(Y,0);
		set_bounds_x(X,shifty); set_bounds_x(Y,shifty);
	} else  {
		if (fjc==1) {
			for (int y=1; y<MY+1; y++) {
				X[0        +y]= Y[BX1*JX+(y+shifty)];
				X[(MX+1)*JX+y]= Y[BXM*JX+(y-shifty)];
				Y[0        +y]= X[BX1*JX+(y+shifty)];
				Y[(MX+1)*JX+y]= X[BXM*JX+(y-shifty)];
			}
		} else {
			cout <<"set_bounds_x error" << endl;
			for (int y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
				for (int k=0; k<fjc; k++) {
					X[k*JX+y]=Y[B_X1[k]*JX+(y+shifty)];
					X[(MX+fjc+k)*JX+y]=Y[B_XM[k]*JX+(y-shifty)];
					Y[k*JX+y]=X[B_X1[k]*JX+(y+shifty)];
					Y[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+(y-shifty)];
				}
			}
		}
	}
}

void LGrad2::set_bounds_y(Real* X, Real*Y, int shiftx){
NAMICS_DBG("set_bounds_y XY in LGrad2 " << endl);
	if (BY1>BYM) {
		//set_bounds_y(X,0); set_bounds_y(Y,0);
		set_bounds_y(X,shiftx); set_bounds_y(Y,shiftx);
	} else {
		if (fjc==1) {
			for (int x=1; x<MX+1; x++) {
				X[x*JX+0   ] =Y[(x+shiftx)*JX+BY1];
				X[x*JX+MY+1] =Y[(x-shiftx)*JX+BYM];
				Y[x*JX+0   ] =X[(x+shiftx)*JX+BY1];
				Y[x*JX+MY+1] =X[(x-shiftx)*JX+BYM];
			}
		} else {
			for (int x=fjc; x<MX+fjc; x++) {
				for (int k=0; k<fjc; k++) {
					X[x*JX+k]=Y[(x+shiftx)*JX+B_Y1[k]];
					X[x*JX+MY+fjc+k]=Y[(x-shiftx)*JX+B_YM[k]];
					Y[x*JX+k]=X[(x+shiftx)*JX+B_Y1[k]];
					Y[x*JX+MY+fjc+k]=X[(x-shiftx)*JX+B_YM[k]];
				}
			}
		}
	}
}

void LGrad2::set_bounds_x(Real* X,int shifty){
NAMICS_DBG("set_bounds_x X in LGrad2 " << endl);	int y;
	int k=0;

	if (fjc==1) {
		for (y=1; y<MY+1; y++) { //not yet 0 and MY+1....
			X[0        +y] = X[BX1*JX+(y+shifty)];
			X[(MX+1)*JX+y] = X[BXM*JX+(y-shifty)];
		}
	} else {
		for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[B_X1[k]*JX+(y+shifty)];
				X[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+(y-shifty)];
			}
		}
	}
}

void LGrad2::set_bounds_y(Real* X,int shiftx){
NAMICS_DBG("set_bounds_y X in LGrad2 " << endl);	int x;
	int k=0;

	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0   ]= X[(x+shiftx)*JX+BY1];
			X[x*JX+MY+1]= X[(x-shiftx)*JX+BYM];
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[(x+shiftx)*JX+B_Y1[k]];
				X[x*JX+MY+fjc+k]=X[(x-shiftx)*JX+B_YM[k]];
			}
		}
	}
}

void LGrad2::set_bounds(Real* X){
NAMICS_DBG("set_bounds in LGrad2 " << endl);	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = X[x*JX+BY1];
			X[x*JX+MY+1]=X[x*JX+BYM];
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = X[BX1*JX+y];
			X[(MX+1)*JX+y]=X[BXM*JX+y];
		}
		//corners....
			X[0] = X[BY1];
			X[MY+1]=X[BYM];
		//}
			X[(MX+1)*JX+      0] = X[(MX+1)*JX+BY1];
			X[(MX+1)*JX+   MY+1]  =X[(MX+1)*JX+BYM];
		//}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+           k]=X[x*JX+B_Y1[k]];
				X[x*JX+    MY+fjc+k]=X[x*JX+B_YM[k]];
			}
		}
		for (y=0; y<MY+2*fjc; y++) {
			for (k=0; k<fjc; k++) {
				X[k*JX+           y]=X[B_X1[k]*JX+     y];
				X[(MX+fjc+k)*JX  +y]=X[B_XM[k]*JX+     y];
			}
		}
			//X[k*JX         +m      ]=X[k*JX            +B_Y1[m]];
			//X[(MX+fjc+k)*JX+m      ]=X[(MX+fjc+k)*JX   +B_Y1[m]];
			//X[k*JX         +(MY+fjc+m)]=X[k*JX         +B_YM[m]];
			//X[(MX+fjc+k)*JX+(MY+fjc+m)]=X[(MX+fjc+k)*JX+B_YM[m]];

			//X[k*JX         +m      ]=X[B_X1[k]*JX+B_Y1[m]];
			//X[(MX+fjc+k)*JX+m      ]=X[B_XM[k]*JX+B_Y1[m]];
			//X[k*JX         +(MY+fjc+m)]=X[B_X1[k]*JX+B_YM[m]];
			//X[(MX+fjc+k)*JX+(MY+fjc+m)]=X[B_XM[k]*JX+B_YM[m]];
		//}
	}
}

void LGrad2::set_M_bounds(Real* X){
NAMICS_DBG("set_bounds in LGrad2 " << endl);	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = X[1*JX+y];
			X[(MX+1)*JX+y]=X[MX*JX+y];
		}
		//corners
		for (x=0; x<1; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
		for (x=MX+1; x<MX+2; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[x*JX+2*fjc-1-k];
				X[x*JX+MY+fjc+k]=X[x*JX+MY+fjc-k-1];
			}
		}
		for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[(2*fjc-1-k)*JX+y];
				X[(MX+fjc+k)*JX+y]=X[(MX+fjc-k-1)*JX+y];
			}
		}
	}
}


void LGrad2::remove_bounds(int *X){
NAMICS_DBG("remove_bounds in LGrad2 " << endl);	int x,y;
	int k;
	if (fjc==1) {
		for (x=0; x<MX+2; x++) {
			X[P(x,0)] = 0;  //needs testing if this is okay to put to zero in case of surface...
			X[P(x,MY+1)]=0;
		}
		for (y=0; y<MY+2; y++) {
			X[P(0,y)] = 0;
			X[P(MX+1,y)]=0;
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=0;
				X[x*JX+MY+fjc+k]=0;
			}
		}
		for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=0;
				X[(MX+fjc+k)*JX+y]=0;
			}
		}
	}
}

void LGrad2::set_bounds(int* X){
NAMICS_DBG("set_bounds in LGrad2 " << endl);	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = X[x*JX+BY1];
			X[x*JX+MY+1]=X[x*JX+BYM];
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = X[BX1*JX+y];
			X[(MX+1)*JX+y]=X[BXM*JX+y];
		}
		//corners
			X[0] = X[BY1];
			X[MY+1]=X[BYM];
		//}
			X[(MX+1)*JX+0] = X[(MX+1)*JX+BY1];
			X[(MX+1)*JX+MY+1]=X[(MX+1)*JX+BYM];
		//}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[x*JX+B_Y1[k]];
				X[x*JX+MY+fjc+k]=X[x*JX+B_YM[k]];
			}
		}
		for (y=fjc; y<MY+fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[B_X1[k]*JX+y];
				X[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+y];
			}
		}
		for (int k=0; k<fjc; k++) for (int m=0; m<fjc; m++)  {
			X[k*JX         +m      ]   =X[B_X1[k]*JX+B_Y1[m]];
			X[(MX+fjc+k)*JX+m      ]   =X[B_XM[k]*JX+B_Y1[m]];
			X[k*JX         +(MY+fjc+m)]=X[B_X1[k]*JX+B_YM[m]];
			X[(MX+fjc+k)*JX+(MY+fjc+m)]=X[B_XM[k]*JX+B_YM[m]];
		}
	}
}

Real LGrad2::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	if (Markov==2) {
		if (lattice_type == hexagonal && !stencil_full) {
			for (int k=0; k<12; k++) {
				GN += WeightedSum(G+k*M);
			}
			GN /=12.0;
			//		GN += 2*WeightedSum(G+k*M);
			//	else
			//		GN += WeightedSum(G+k*M);
			//}
			//GN /=12.0;

		} else {
			if (lattice_type==simple_cubic) {
				for (int k=0; k<5; k++) {
					if (k== 2)
						GN += 2*WeightedSum(G+k*M);
					else
						GN += WeightedSum(G+k*M);
				}
				GN /=6.0;
			} else {
				for (int k=0; k<12; k++) {
					GN += WeightedSum(G+k*M);
				}
				GN /=12.0;
			}
		}
	} else GN = WeightedSum(G);
	return GN;
}
void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb,int Markov, int M){
	if (Markov==2) {
		if (lattice_type == hexagonal&& !stencil_full) {
			Real C=1.0/12.0;
			for (int k=0; k<12; k++) {
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
			}
			//	else
			//}

		} else {
			if (lattice_type==simple_cubic) {
				Real C1=1.0/3.0;
				Real C2=1.0/6.0;
				for (int k=0; k<5; k++) {
					if (k==2)
						for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C1) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
					else
						for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C2) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
				}
			} else { //hexagonal
				Real C=1.0/12.0;
				for (int k=0; k<12; k++) {
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (C) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
				}
			}
		}
	} else for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (Gf)[__i] * (Gb)[__i];
}
void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb,Real degeneracy, int Markov, int M){
	if (Markov==2) {
		if (lattice_type == hexagonal&& !stencil_full) {
			for (int k=0; k<12; k++) {
				for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/12.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
			}
			//	else {
			//	}
			//}

		} else {
			if (lattice_type==simple_cubic) {
				for (int k=0; k<5; k++)
				if (k==2)
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/3.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
				else
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/6.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
			} else {
				for (int k=0; k<12; k++) {
					for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/12.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
				}
			}
		}
	} else for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy) * (Gf)[__i] * (Gb)[__i];
}


void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb, Real* G1, Real norm, int Markov, int M){
	(void)M;
	(void)Markov;
	(void)norm;
	(void)G1;
	(void)Gb;
	(void)Gf;
	(void)phi;
	cout << "composition not (yet) implemented for alias in LGrad2 " << endl;
}

void LGrad2::Initiate(Real* G,Real* Gz,int Markov, int M){
	if (Markov==2) {
		if (lattice_type == hexagonal&& !stencil_full) {
			for (int k=0; k<12; k++) std::copy_n(Gz, M, G+k*M);
		} else {
			if (lattice_type==simple_cubic) {
				for (int k=0; k<5; k++) std::copy_n(Gz, M, G+k*M);
			} else {
				for (int k=0; k<12; k++) std::copy_n(Gz, M, G+k*M);
			}
		}
	} else std::copy_n(Gz, M, G);
}

void LGrad2::Terminate(Real* Gz,Real* G,int Markov, int M){
NAMICS_DBG("LGrad2:: terminate " << endl);	if (Markov==2) {
		cout <<"terminate in markov==2 is not tested" << endl;
		std::fill_n(Gz, M, 0);
		if (lattice_type == hexagonal&& !stencil_full) {
			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (Gz)[__i] += (G+k*M)[__i];
		} else {
			Real C=1.0/6.0;
			if (lattice_type==simple_cubic) {
				for (int k=0; k<5; k++) for (int __i = 0; __i < (M); ++__i) (Gz)[__i] += (G+k*M)[__i];
				for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (C);  //dit is in ieder geval ook niet goed.
			} else {
				for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (Gz)[__i] += (G+k*M)[__i];
				for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (C);
			}
		}
		cout <<"possible problem in LGrad2::Terminate " << endl;
	} else std::copy_n(G, M, Gz);
}

bool LGrad2:: PutMask(Real* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	(void)pz;
NAMICS_DBG("PutMask in LGrad2 " << endl);	//R*=fjc; //is already done in segment
	bool success=true;
	int length =px.size();
	int X,Y;
	int dx,dy;
	Real teller,noemer;
	if (length > 1) {
		cout <<"In two gradient system, we can have just one particle: we found " <<length <<"particles. " << endl;
		return false;
	}
	for (int i =0; i<length; i++) {
		int xx,yy;
		xx=px[i]; yy=py[i];
		if (xx !=0) {
			cout <<"In two gradients system, we expect the particle at the central axis" << endl;
			return false;
		}
		if (R>MX || R>MY) {cout <<" particle should be smaller than size of box in X or Y direction" << endl; return false;}
		for (int x=1; x<R+2; x++)
		for (int y=yy-R; y<yy+R+2; y++){
			X=x; Y=y;
			if (x*x+(yy-y)*(yy-y) <=(R+1)*(R+1)) {
				if (x*x+(yy-y)*(yy-y) <=(R-1)*(R-1)) {
					if (!(y<fjc || y>MY+fjc-1))  MASK[P(X,Y)]++;
				} else {
					teller=0; noemer=0;
					for (dx=0; dx<10; dx++) {
						for (dy=0; dy<10; dy++) {
							noemer +=x+dx/10;
							if ((x+dx)*(x+dx)+(yy-y-dy)*(yy-y-dy) <=R*R) teller +=x+dx/10;
						}
					}
					if (!(y<fjc || y>MY+fjc-1))  MASK[P(X,Y)]=teller/noemer;
				} //at the edge
			} //else too large
		}
	}
	return success;
}


Real LGrad2::DphiDt(Real* g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB,Real B_A, Real B_B) {
	NAMICS_DBG("LGrad2: DphiDt not implemented yet " << endl);	Real AverageJ=0;
	int x, y;

	//apply BC to the edges
	for (x=1; x<MX+1; x++) {
		//lower edge = south/current-1
		g[x*JX+1] = phiA[x*JX+0]/phiA[x*JX+1] - 1.0;
		//upper edge = north/current-1
		g[x*JX+MY] = phiA[x*JX+(MY+1)]/phiA[x*JX+MY] - 1.0;
	}
	for (y=1; y<MY+1; y++) {
		//left edge = west/current-1
		g[1*JX + y] =  phiA[0*JX+y]/phiA[1*JX+y] - 1.0;
		//right edge = east/current -1
		g[MX*JX + y] =  phiA[(MX+1)*JX+y]/phiA[MX*JX+y] - 1.0;
	}

	//gets indices for the neighbors
	int north, south, east, west, center;
	auto get_c = [&](){return x*JX +y;};
	auto get_w = [&](){return (x-1)*JX +y;};
	auto get_e = [&](){return (x+1)*JX +y;};
	auto get_s = [&](){return x*JX +(y-1);};
	auto get_n = [&](){return x*JX +(y+1);};
	auto update_neighbors = [&](){
		center = get_c();
		north = get_n();
		south = get_s();
		east = get_e();
		west = get_w();
	};

	//memory efficient, but computationally inefficient
	//the flux has to be calculated first
	//and only north and east component stored
	Real a,b,c,Ma,Mb,Mc;
	for (y=2; y<MY; y++){
		for (x=2; x<MX; x++){
			update_neighbors();

			a = phiA[west]*phiB[west]*B_B/B_phitot[west];
			b = phiA[center]*phiB[center]*B_B/B_phitot[center];
			c = phiA[east]*phiB[east]*B_B/B_phitot[east];

			Ma = alphaA[west]-alphaB[west];
			Mb = alphaA[center]-alphaB[center];
			Mc = alphaA[east]-alphaB[east];

			//flux divergence x-axis
			g[x*JX + y] += (a+b)*(Mb - Ma)*lambda_1[center] - (b+c)*(Mc-Mb)*lambda1[center];
			AverageJ += (a+b)*(Mb - Ma)*lambda_1[center]*L[center];//L[x]?

			a = phiA[south]*phiB[south]*B_B/B_phitot[south];
			c = phiA[north]*phiB[north]*B_B/B_phitot[north];

			Ma = alphaA[south]-alphaB[south];
			Mc = alphaA[north]-alphaB[north];

			//flux_divergence y-axis
			g[x*JX + y] += (a+b)*(Mb - Ma)*lambda0[center] - (b+c)*(Mc-Mb)*lambda0[center];
			AverageJ+=(b+c)*(Mc-Mb)*lambda0[center];
		}
	}


	return -B_A*AverageJ/(2*(M-MX*2-MY*2)*lambda);

}

Real LGrad2::MomentPlanar(Real* X,int n,Real Z0){
	(void)Z0;
	(void)n;
	(void)X;
	cout <<"MomentPlanar not implemented; kJ0 or kbar may be wrong. " << endl;
	return 0;
}
