#include <iostream>
#include <string>
#include "lattice.h"
#include "LG2Planar.h"

LG2Planar::LG2Planar(const Input& In_,const string& name_): LGrad2(In_,name_) {
	JY=1;
}

LG2Planar::~LG2Planar() {
NAMICS_DBG("LG2Planar destructor " << endl);}


void LG2Planar:: ComputeLambdas() {

	if (fjc==1) {
		for (int i=0; i<M; i++) L[i]=1;
	}

	if (fjc==2) {
		for (int i=0; i<M; i++) L[i]=1.0/fjc;
	}
}


void LG2Planar::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG(" Side in LG2Planar " << endl);	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0);set_bounds(X);

	if (fcc_sites) {
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X+1)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+1)[__i] += (X+JX)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JX+1)[__i];
		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (1.0/9.0);
	} else {
		if (fjc==1) {
			if (!stencil_full) {
				if (lattice_type==simple_cubic) { //6 point stencil (voorheen 9-punts
					for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
					for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i];
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (1.0/2.0);
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i];
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (1.0/3.0);
				} else { //Johan's method
					for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
					for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (X)[__i];
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (2.0);

					remove_bounds(X);
					set_bounds_x(X,-1);
					for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JX)[__i] += (X+JY)[__i];
					for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JY)[__i] += (X+JX)[__i];

					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (1.0/12.0);

					}
			} else {
				if (lattice_type==simple_cubic) {  //9punts stencil, voorheen 6pnts
					Real C1=16.0/36.0;
					Real C2=4.0/36.0;
					Real C3=1.0/36.0;
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (C1) * (X)[__i];
					for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (C2) * (X)[__i];
					for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (C2) * (X+1)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (C2) * (X)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (C2) * (X+JX)[__i];
					for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (C3) * (X)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (C3) * (X+1)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side+1)[__i] += (C3) * (X+JX)[__i];
					for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (C3) * (X+JX+1)[__i];

				} else {
						//hexagonal //9 point stencil
					Real Two=2.0;
					Real C=1.0/16.0;
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (Two) * (X)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
					for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Two);
					for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JX+JY)[__i] += (X)[__i];
					for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side)[__i] += (X+JX+JY)[__i];
					for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JX)[__i] += (X+JY)[__i];
					for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JY)[__i] += (X+JX)[__i];

					for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);
				}
			}

		} else {  //fjc>1
			for (int block=0; block<3; block++){
				int bk;
				int a,b;
				for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) {
					bk=a=b=0;
					if (x==-fjc || x==fjc) bk++;
					if (y==-fjc || y==fjc) bk++;
					if (bk==block) {
						if (x<0) a =-x*JX; else b=x*JX;
						if (y<0) a -=y*JY; else b+=y*JY;
						for (int __i = 0; __i < (M-a-b); ++__i) (X_side+a)[__i] += (X+b)[__i];
					}
				}
				if (block !=2) for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (2.0); else for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (1.0/(4.0*(FJC-2)*FJC+1));
			}
		}
	}
}

void LG2Planar::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {
		if (lattice_type==hexagonal) {

			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs,   *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0);set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);
			set_bounds_y(gz2,gz9,0);  set_bounds_y(gz3,gz8,0); set_bounds_y(gz4,gz7,0); set_bounds_y(gz0,gz11,0); set_bounds_y(gz1,gz10,0); set_bounds_y(gz5,gz6,0);

			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz0)[__i]; //0 and 1 are equivalent
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz1)[__i]; //10 and 11 are equivalent
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz2)[__i]; //3 and 4 equivalent
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz3)[__i]; //7 and 8
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz8)[__i];


			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz3+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz4+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz5+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz6+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz7+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[1]) * (gz8+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[0]) * (gz9+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[0]) * (gz10+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11)[__i] += (P[0]) * (gz11+JX)[__i];


			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_y(gz2,gz9,0);  set_bounds_y(gz3,gz8,0); set_bounds_y(gz4,gz7,0); set_bounds_y(gz0,gz11,0); set_bounds_y(gz1,gz10,0); set_bounds_y(gz5,gz6,0);

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

			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

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

			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[0]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz8)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz3+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz4+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz5+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz6+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz7+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[1]) * (gz8+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[0]) * (gz9+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[0]) * (gz10+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10)[__i] += (P[0]) * (gz11+JX)[__i];

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

			set_bounds_y(gz0,gz11,0);
			set_bounds_y(gz1,gz10,0);
			set_bounds_y(gz2,gz9,0);
			set_bounds_y(gz3,gz8,0);
			set_bounds_y(gz4,gz7,0);
			set_bounds_y(gz5,gz6,0);

			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz5+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz6+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz7+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz8+JY)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz9)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz11)[__i];

			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];

		} else { //simple_cubic
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
			Real *g=G1;

			std::fill_n(gs, 5*M, 0);
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz0)[__i];
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
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[0]) * (gz4+JX)[__i];

			for (int k=0; k<5; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		if (lattice_type == simple_cubic) { //size=2*FJC-1 here FJC=3
			if (fjc==1) {

			} else { //fjc==2

			/* geheugensteun
			0   x   x
			1   x   y
			2   x  -y
			3   x   0
			4   y   x
			5   y   y
			6   y   -x
			7   y  0
			8   0   y
			9   0   x
			10  0   0
			11  0  -x
			12  0  -y
			13 -y   0
			14 -y  x
			15 -y  -y
			16 -y -x
			17 -x   0
			18 -x  -y
			19 -x   y
			20 -x  -x
			*/
				int size = 2*FJC-1;
				Real *gs=G+M*size*s_to;
				Real *gs_1=G+M*size*s_from;
				Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
				Real            *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M;
				set_bounds_x(gz0,gz8,0);
				set_bounds_x(gz1,gz7,0);
				set_bounds_x(gz2,gz6,0);
				set_bounds_x(gz3,0);
				set_bounds_x(gz4,0);
				set_bounds_x(gz5,0);

				set_bounds_y(gz1,gz6,0);
				set_bounds_y(gz2,gz7,0);
				set_bounds_y(gz3,gz5,0);
				set_bounds_y(gz0,0);
				set_bounds_y(gz5,0);
				set_bounds_y(gz8,0);
				Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
				Real *g=G1;

				std::fill_n(gs, size*M, 0);
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[0]) * (gz0)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[1]) * (gz1)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[1]) * (gz2)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[2]) * (gz3)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[2]) * (gz4)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[2]) * (gz5)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[3]) * (gz6)[__i];
				for (int __i = 0; __i < (M-2*JX); ++__i) (gx0+2*JX)[__i] += (P[3]) * (gz7)[__i];

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
				for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[0]) * (gz4+JX)[__i];
				for (int k=0; k<size; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
			}
		} else {//hexagonal // not correct because molecules can not 'diffuse' from one plane to the other.... ('behoud' wel segmenten....)
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M;
			Real *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M;
			Real *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			set_bounds_x(gz0,gz3,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0);
			set_bounds_y(gz1,gz2,0); set_bounds_y(gz0,0); set_bounds_y(gz3,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M;
			Real *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M;
			Real *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz2)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz3)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[0]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz3+JY)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[0]) * (gz3+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[0]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (gz6)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[1]) * (gz7)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[1]) * (gz7)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[1]) * (gz5+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[1]) * (gz6+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[0]) * (gz7+JX)[__i];


			for (int __i = 0; __i < (M); ++__i) (gx8)[__i] += (P[0]) * (gz8)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx8)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx8)[__i] += (P[1]) * (gz10)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz8)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz9)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz11)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[1]) * (gz8+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[0]) * (gz10+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[1]) * (gz11+JY)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx11)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx11)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx11)[__i] += (P[0]) * (gz11)[__i];

			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];		}
	}
}
void LG2Planar::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {
		if (lattice_type==hexagonal) { //0-2 9-11


			Real *gs=G+M*12*s_to;  //3-4
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			for (int k=0; k<12; k++) remove_bounds(gs_1+k*M);

			set_bounds_x(gz0,gz11,0);

			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5+JX)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6+JX)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7+JX)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8+JX)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx9+JX)[__i] += (P[0]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10+JX)[__i] += (P[0]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11+JX)[__i] += (P[0]) * (gz11)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[0]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[0]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8)[__i] += (P[1]) * (gz0+JX)[__i];

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

			set_bounds_x(gz1,gz10,0);

			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5+JX)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6+JX)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7+JX)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8+JX)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx9+JX)[__i] += (P[0]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx10+JX)[__i] += (P[0]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx11+JX)[__i] += (P[0]) * (gz10)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[0]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[0]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx7)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx8)[__i] += (P[1]) * (gz1+JX)[__i];

			remove_bounds(gz1);remove_bounds(gz10);
			set_bounds_y(gz2,gz9,0);

			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx6)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9)[__i] += (P[0]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10)[__i] += (P[0]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11)[__i] += (P[0]) * (gz9+JY)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx6+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx7+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8+JY)[__i] += (P[1]) * (gz2)[__i];

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

		} else {
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M,*gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
			Real *g=G1;

			std::fill_n(gs, 5*M, 0);
			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2+JX)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[0]) * (gz4)[__i];

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

			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (gz0+JX)[__i];

			for (int k=0; k<5; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];		}

	} else {
		if (lattice_type==simple_cubic) {
			cout <<"simple_cubic, markov 2, planar, stencil_full, not implemented " << endl;
		} else { //hexagonal //not finished and also not physically acceptable....there are three sub-lattices and there is no cross-over from one to the other.
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M;
			Real *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M;
			Real *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			set_bounds_x(gz0,gz3,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0);
			set_bounds_y(gz1,gz2,0); set_bounds_y(gz0,0); set_bounds_y(gz3,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M;
			Real *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M;
			Real *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx0)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (gz0+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1)[__i] += (P[0]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[1]) * (gz0+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx2+JX)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[1]) * (gz0+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[1]) * (gz1+JY)[__i];


			for (int __i = 0; __i < (M); ++__i) (gx4)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx4)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[0]) * (gz4+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx5+JX)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx5)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz4+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx6+JX)[__i] += (P[1]) * (gz7)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx6)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx6)[__i] += (P[1]) * (gz4+JX)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx7+JX)[__i] += (P[0]) * (gz7)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx7)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx7)[__i] += (P[1]) * (gz5)[__i];


			for (int __i = 0; __i < (M-JY); ++__i) (gx8+JY)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx8)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx8)[__i] += (P[0]) * (gz8)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx9)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx9)[__i] += (P[0]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx9)[__i] += (P[1]) * (gz8)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx10)[__i] += (P[1]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx10+JY)[__i] += (P[0]) * (gz10)[__i];
			for (int __i = 0; __i < (M); ++__i) (gx10)[__i] += (P[1]) * (gz8)[__i];

			for (int __i = 0; __i < (M); ++__i) (gx11)[__i] += (P[0]) * (gz11)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11+JY)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx11)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];		}
	}
}


void LG2Planar::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
NAMICS_DBG(" propagate in LGrad2 " << endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	std::fill_n(gs, M, 0); set_bounds(gs_1);
	if (fjc==1) {
		if (!stencil_full) {
			if (lattice_type==simple_cubic) { //9 point stencil..
				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/2.0);
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/3.0);
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];


			} else { //hexagonal Johan's method //kept for nostalgic reasons

				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (2.0);

				remove_bounds(gs_1);
				set_bounds_x(gs_1,-1);
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JX)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JY)[__i] += (gs_1+JX)[__i];

				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/12.0);
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
				}
		} else {
			if (lattice_type==simple_cubic) {
				Real C1=16.0/36.0;
				Real C2=4.0/36.0;
				Real C3=1.0/36.0;
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (C1) * (gs_1)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (C2) * (gs_1)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (C2) * (gs_1+1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (C2) * (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (C2) * (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX+1)[__i] += (C3) * (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (C3) * (gs_1+1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs+1)[__i] += (C3) * (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (gs)[__i] += (C3) * (gs_1+JX+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
			} else { //hexagonal //9 point stencil

				Real Two=2.0;
				Real C=1.0/16.0;
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (Two) * (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (Two);
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JX+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs)[__i] += (gs_1+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JX)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JY)[__i] += (gs_1+JX)[__i];

				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C);
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
			}
		}
	} else {
		for (int block=0; block<3; block++){
			Real Two=2.0;
			int bk;
			int a,b;
			for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) {
				bk=a=b=0;
				if (x==-fjc || x==fjc) bk++;
				if (y==-fjc || y==fjc) bk++;
				if (bk==block) {
					if (x<0) a =-x*JX; else b=x*JX;
					if (y<0) a -=y*JY; else b+=y*JY;
					for (int __i = 0; __i < (M-a-b); ++__i) (gs+a)[__i] += (gs_1+b)[__i];
				}
			}
			if (block !=2) for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (Two); else for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/(4.0*(FJC-2)*FJC+1));
		}
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
	}
}


void LG2Planar::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	std::fill_n(EE, M, 0);
	Real Exmin,Explus,Eymin,Eyplus;
	int x,y,z;

	pf = pf/2.0;
	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++) {
			z=x*JX+y;
			Exmin=psi[z]-psi[z-JX];
			Exmin*=Exmin;
			Explus=psi[z]-psi[z+JX];
			Explus*=Explus;
			Eymin=psi[z]-psi[z-1];
			Eymin*=Eymin;
			Eyplus=psi[z]-psi[z+1];
			Eyplus*=Eyplus;
			EE[x*JX+y]=pf*(Exmin+Explus+Eymin+Eyplus);
		}
	}
}


void LG2Planar::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, Real* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y,i;

	Real epsXplus, epsXmin, epsYplus,epsYmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);
	Real ax,ay;

	if (!fixedPsi0) {
		C=C*2.0/fjc/fjc;
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				i=x*JX+y;
				epsXmin=eps[i]+eps[i-JX];
				epsXplus=eps[i]+eps[i+JX];
				epsYmin=eps[i]+eps[i-1];
				epsYplus=eps[i]+eps[i+1];
				if (x==fjc) ax=psi[i-JX]; else ax=X[i-JX]; //upwind
				if (y==fjc) ay=psi[i-1]; else ay=X[i-1]; //upwind
				//X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/(epsXmin+epsXplus+epsYmin+epsYplus);
				X[i]= (C*q[i]+epsXmin*ax+epsXplus*psi[i+JX]+epsYmin*ay+epsYplus*psi[i+1])/(epsXmin+epsXplus+epsYmin+epsYplus);
			}
		}
		for (int __i = 0; __i < (M); ++__i) (g)[__i] = (g)[__i] - (X)[__i];
  	} else { //fixedPsi0 is true
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++){
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]=0.25*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
			        +psi[x*JX+y-1]  +psi[x*JX+y+1])
				 +0.5*q[x*JX+y]*C/eps[x*JX+y];
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


void LG2Planar::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, Real* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
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

bool LG2Planar:: PutMask(Real* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	(void)R;
	(void)pz;
	(void)py;
	(void)px;
	(void)MASK;
	bool success=false;
	cout <<"PutMask does not make sence in planar 2 gradient system " << endl;
	return success;
}
