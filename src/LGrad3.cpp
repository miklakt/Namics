#include <iostream>
#include <string>
#include "lattice.h"
#include "LGrad3.h"

LGrad3::LGrad3(const Input& In_,const string& name_): Lattice(In_,name_) {}

LGrad3::~LGrad3() {
NAMICS_DBG_THIS("LGrad3 destructor " << endl);}

void LGrad3:: ComputeLambdas() {
}

bool LGrad3::PutM() {
NAMICS_DBG_THIS("PutM in LGrad3 " << endl);	bool success=true;
	volume = MX*MY*MZ;
	JX=(MZ+2*fjc)*(MY+2*fjc); JY=MZ+2*fjc; JZ=1; M = (MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc);

	Accesible_volume=volume;
	return success;
}

void LGrad3::TimesL(Real* X){
	(void)X;
NAMICS_DBG_THIS("TimesL in LGrad3 " << endl);}

void LGrad3::DivL(Real* X){
	(void)X;
NAMICS_DBG_THIS("DivL in LGrad3 " << endl);}

Real LGrad3:: Moment(Real* X,Real Xb, int n) {
	(void)n;
	(void)Xb;
	(void)X;
NAMICS_DBG_THIS("Moment in LGrad3 " << endl);	Real Result=0;
	return Result/fjc;
}

Real LGrad3::WeightedSum(Real* X){
NAMICS_DBG_THIS("weighted sum in LGrad3 " << endl);	Real sum{0};
	remove_bounds(X);
	(sum) = 0; for (int __i = 0; __i < (M); ++__i) (sum) += (X)[__i];
	return sum;
}

void LGrad3::vtk(string filename, Real* X, string id,bool writebounds) {
	(void)filename;
	(void)X;
	(void)id;
	(void)writebounds;
NAMICS_DBG_THIS("vtk in LGrad3 " << endl);	cout << "VTK output is disabled; use kal/pro output instead." << endl;
}

void LGrad3::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds,bool DOS){
NAMICS_DBG_THIS("PutProfiles in LGrad3 " << endl);	Real one=1.0;
	int x,y,z,i;
	int length=X.size();
	int a;
	if (writebounds) a=0; else a = fjc;
	for (x=a; x<MX+2*fjc-a; x++)
	for (y=a; y<MY+2*fjc-a; y++)
	for (z=a; z<MZ+2*fjc-a; z++) {
#ifdef LongReal
		writer->Writef(pf,"%Le\t%Le\t%Le\t",one*(x-fjc+1)/fjc-0.5/fjc,one*(y-fjc+1)/fjc-0.5/fjc,one*(z-fjc+1)/fjc-0.5/fjc);
		for (i=0; i<length; i++) writer->Writef(pf,"%.20Le\t",X[i][P(x,y,z)]);
#else
		writer->Writef(pf,"%e\t%e\t%e\t",one*(x-fjc+1)/fjc-0.5/fjc,one*(y-fjc+1)/fjc-0.5/fjc,one*(z-fjc+1)/fjc-0.5/fjc);
		for (i=0; i<length; i++) writer->Writef(pf,"%.20e\t",X[i][P(x,y,z)]);

#endif
		if (DOS) writer->Writef(pf,"\r\n"); else writer->Writef(pf,"\n");
	}
}

void LGrad3::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG_THIS(" Side in LGrad3 " << endl);	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0);//set_bounds(X);
	Real Two=2.0;
	Real Four=4.0;
	if (!stencil_full) {
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
		for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i];

		if (lattice_type == simple_cubic) {
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Four);
		} else {
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Two);
		}
		for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JX+JY)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side)[__i] += (X+JX+JY)[__i];
		for (int __i = 0; __i < (M-JY-JX); ++__i) (X_side+JY)[__i] += (X+JX)[__i];
		for (int __i = 0; __i < (M-JY-JX); ++__i) (X_side+JX)[__i] += (X+JY)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JX+1)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X+1)[__i];
		for (int __i = 0; __i < (M-JX); ++__i) (X_side+1)[__i] += (X+JX)[__i];
		for (int __i = 0; __i < (M-JY-1); ++__i) (X_side+JY+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JY+1)[__i];
		for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X+1)[__i];
		for (int __i = 0; __i < (M-JY); ++__i) (X_side+1)[__i] += (X+JY)[__i];
		if (lattice_type == simple_cubic) {
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Four);
		} else {
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Two);
		}
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+JY+1)[__i] += (X)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side)[__i] += (X+JX+JY+1)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+JY)[__i] += (X+1)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+1)[__i] += (X+JX+JY)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+1)[__i] += (X+JY)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JY)[__i] += (X+JX+1)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JY+1)[__i] += (X+JX)[__i];
		for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX)[__i] += (X+JY+1)[__i];
		if (lattice_type == simple_cubic) {
			Real C=1.0/152.0;
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);
		} else {
			Real C=1.0/56.0;
			for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);
		}
	} else {
		if (lattice_type==simple_cubic) {
			Real C=1.0/6.0;
			for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (X_side+JZ)[__i] += (X)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (X_side)[__i] += (X+JZ)[__i];
	 		for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);
		} else { //hexagonal
			if (fjc==1) {
				Real Two=2.0;
				Real C=1.0/40.0;
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] += (Two) * (X)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Two);

				for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JX+JY)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side)[__i] += (X+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JX)[__i] += (X+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (X_side+JY)[__i] += (X+JX)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX+1)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+1)[__i] += (X+JX)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (X_side+JY+1)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (X_side+1)[__i] += (X+JY)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side+JX)[__i] += (X+1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (X_side)[__i] += (X+JX+1)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (X_side+JY)[__i] += (X+1)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (X_side)[__i] += (X+JY+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Two);

				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+JY+1)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side)[__i] += (X+JX+JY+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+JY)[__i] += (X+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+1)[__i] += (X+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+1)[__i] += (X+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JY)[__i] += (X+JX+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JY+1)[__i] += (X+JX)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX)[__i] += (X+JY+1)[__i];

				for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);


			} else { //fjc==2
				for (int block=0; block<4; block++){
					Real Two=2.0;
					Real C=1.0/8.0/((FJC-2)*(FJC-2)*(FJC-2)+3*(FJC-2)*(FJC-2)+3*(FJC-2)+1);
					int bk;
					int a,b;
					for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) for (int z=-fjc; z<fjc+1; z++){
						bk=a=b=0;
						if (x==-fjc || x==fjc) bk++;
						if (y==-fjc || y==fjc) bk++;
						if (z==-fjc || z==fjc) bk++;
						if (bk==block) {
							if (x<0) a =-x*JX; else b=x*JX;
							if (y<0) a -=y*JY; else b+=y*JY;
							if (z<0) a -=z*JZ; else b+=z*JZ;
							for (int __i = 0; __i < (M-a-b); ++__i) (X_side+a)[__i] += (X+b)[__i];
						}
					}
					if (block !=3) for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (Two); else for (int __i = 0; __i < (M); ++__i) (X_side)[__i] *= (C);

				}
			}
		}
	}
}

void LGrad3::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {

		if (lattice_type == hexagonal ) {


			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs,   *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0,0); set_bounds_x(gz1,gz10,0,0);set_bounds_x(gz2,gz9,0,0); set_bounds_x(gz3,gz8,0,0); set_bounds_x(gz4,gz7,0,0); set_bounds_x(gz5,gz6,0,0);

			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz3)[__i];
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

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz2,gz9,0,0);  set_bounds_y(gz3,gz8,0,0); set_bounds_y(gz4,gz7,0,0); set_bounds_y(gz0,gz11,0,0); set_bounds_y(gz1,gz10,0,0); set_bounds_y(gz5,gz6,0,0);

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
			set_bounds_z(gz1,gz10,0,0); set_bounds_z(gz4,gz7,0,0); set_bounds_z(gz5,gz6,0,0); set_bounds_z(gz0,gz11,0,0); set_bounds_z(gz2,gz9,0,0); set_bounds_z(gz3,gz8,0,0);

			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[0]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[0]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz9)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz10)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz11)[__i];

			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[1]) * (gz0+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[1]) * (gz1+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[1]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[0]) * (gz6+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[0]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[0]) * (gz8+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[1]) * (gz9+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6)[__i] += (P[1]) * (gz11+JZ)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0,-1); set_bounds_x(gz1,gz10,0,-1);set_bounds_x(gz2,gz9,0,-1); set_bounds_x(gz3,gz8,0,-1); set_bounds_x(gz4,gz7,0,-1); set_bounds_x(gz5,gz6,0,-1);

			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[0]) * (gz0+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[0]) * (gz1+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[0]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[1]) * (gz3+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[1]) * (gz4+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[1]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[1]) * (gz6+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[1]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx1+JX)[__i] += (P[1]) * (gz8+JZ)[__i];

			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz3+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz4+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz5+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz6+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz7+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz8+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[0]) * (gz9+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[0]) * (gz10+JX)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JZ)[__i] += (P[0]) * (gz11+JX)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz2,gz9,-1,0);  set_bounds_y(gz3,gz8,-1,0); set_bounds_y(gz4,gz7,-1,0); set_bounds_y(gz0,gz11,-1,0); set_bounds_y(gz1,gz10,-1,0); set_bounds_y(gz5,gz6,-1,0);

			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[0]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[0]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[0]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[1]) * (gz5+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[1]) * (gz6+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[1]) * (gz7+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx2+JX)[__i] += (P[1]) * (gz8+JY)[__i];

			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz3+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz4+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz5+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz6+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz7+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[1]) * (gz8+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz9+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz10+JX)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JY)[__i] += (P[0]) * (gz11+JX)[__i];

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_z(gz1,gz10,0,-1); set_bounds_z(gz4,gz7,0,-1); set_bounds_z(gz5,gz6,0,-1); set_bounds_z(gz0,gz11,0,-1); set_bounds_z(gz2,gz9,0,-1); set_bounds_z(gz3,gz8,0,-1);


			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[1]) * (gz0+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[1]) * (gz1+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[1]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[0]) * (gz3+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[0]) * (gz4+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[0]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[1]) * (gz9+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JY)[__i] += (P[1]) * (gz11+JZ)[__i];

			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[0]) * (gz6+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[0]) * (gz7+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[0]) * (gz8+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz10+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz11+JY)[__i];

			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		} else {
			Real *gs=G+M*6*s_to;
			Real *gs_1=G+M*6*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M;
			set_bounds_x(gz0,gz5,0,0); set_bounds_x(gz1,0,0); set_bounds_x(gz2,0,0); set_bounds_x(gz3,0,0); set_bounds_x(gz4,0,0);
			set_bounds_y(gz1,gz4,0,0); set_bounds_y(gz0,0,0); set_bounds_y(gz2,0,0); set_bounds_y(gz3,0,0); set_bounds_y(gz5,0,0);
			set_bounds_z(gz2,gz3,0,0); set_bounds_z(gz0,0,0); set_bounds_z(gz1,0,0); set_bounds_z(gz4,0,0); set_bounds_z(gz5,0,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M;
			Real *g=G1;

			std::fill_n(gs, 6*M, 0);
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[0]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx0+JX)[__i] += (P[1]) * (gz4)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[0]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1+JY)[__i] += (P[1]) * (gz5)[__i];

			for (int __i = 0; __i < (M-JZ); ++__i) (gx2+JZ)[__i] += (P[1]) * (gz0)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2+JZ)[__i] += (P[1]) * (gz1)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2+JZ)[__i] += (P[0]) * (gz2)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2+JZ)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2+JZ)[__i] += (P[1]) * (gz5)[__i];

			for (int __i = 0; __i < (M-JZ); ++__i) (gx3)[__i] += (P[1]) * (gz0+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx3)[__i] += (P[1]) * (gz1+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx3)[__i] += (P[0]) * (gz3+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx3)[__i] += (P[1]) * (gz4+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx3)[__i] += (P[1]) * (gz5+JZ)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (gz0+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (gz2+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (gz3+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4)[__i] += (P[1]) * (gz5+JY)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz3+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[1]) * (gz4+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5)[__i] += (P[0]) * (gz5+JX)[__i];

			for (int k=0; k<6; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		if (lattice_type ==simple_cubic) {
			cout <<"Markov 2, simple cubic, geometry 3, stencil_full not implemented" << endl;
		} else {
			cout <<"Markov 2, hexagonal, geometry 3, stencil_full not implemented" << endl;
		}
	}
}
void LGrad3::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {
		if (lattice_type == hexagonal) {


			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			std::fill_n(gs, 12*M, 0);
			for (int k=0; k<12; k++) remove_bounds(gs_1+k*M);
			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);
			//remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0,0);

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
			set_bounds_y(gz3,gz8,0,0);

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
			set_bounds_z(gz5,gz6,0,0);

			for (int __i = 0; __i < (M-JZ); ++__i) (gx0+JZ)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx1+JZ)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2+JZ)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx6+JZ)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx7+JZ)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx8+JZ)[__i] += (P[0]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx9+JZ)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz6)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx11+JZ)[__i] += (P[1]) * (gz6)[__i];

			for (int __i = 0; __i < (M-JZ); ++__i) (gx0)[__i] += (P[1]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx1)[__i] += (P[1]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2)[__i] += (P[1]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx3)[__i] += (P[0]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx4)[__i] += (P[0]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5)[__i] += (P[0]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx9)[__i] += (P[1]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx10)[__i] += (P[1]) * (gz5+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx11)[__i] += (P[1]) * (gz5+JZ)[__i];

			remove_bounds(gz5); remove_bounds(gz6);
			set_bounds_x(gz1,gz10,0,-1);

			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx3+JX)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx4+JX)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx5+JX)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx6+JX)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx7+JX)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx8+JX)[__i] += (P[1]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx9+JX)[__i] += (P[0]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx10+JX)[__i] += (P[0]) * (gz10+JZ)[__i];
			for (int __i = 0; __i < (M-JX-JZ); ++__i) (gx11+JX)[__i] += (P[0]) * (gz10+JZ)[__i];

			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx0+JZ)[__i] += (P[0]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx1+JZ)[__i] += (P[0]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx2+JZ)[__i] += (P[0]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx3+JZ)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx4+JZ)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx6+JZ)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx7+JZ)[__i] += (P[1]) * (gz1+JX)[__i];
			for (int __i = 0; __i < (M-JZ-JX); ++__i) (gx8+JZ)[__i] += (P[1]) * (gz1+JX)[__i];

			remove_bounds(gz1);remove_bounds(gz10);
			set_bounds_y(gz2,gz9,-1,0);

			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx3+JX)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx4+JX)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx5+JX)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx6+JX)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx7+JX)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx8+JX)[__i] += (P[1]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx9+JX)[__i] += (P[0]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx10+JX)[__i] += (P[0]) * (gz9+JY)[__i];
			for (int __i = 0; __i < (M-JX-JY); ++__i) (gx11+JX)[__i] += (P[0]) * (gz9+JY)[__i];

			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx0+JY)[__i] += (P[0]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx1+JY)[__i] += (P[0]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx2+JY)[__i] += (P[0]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx3+JY)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx4+JY)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx5+JY)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx6+JY)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx7+JY)[__i] += (P[1]) * (gz2+JX)[__i];
			for (int __i = 0; __i < (M-JY-JX); ++__i) (gx8+JY)[__i] += (P[1]) * (gz2+JX)[__i];

			remove_bounds(gz2);remove_bounds(gz9);
			set_bounds_z(gz4,gz7,0,-1);

			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx0+JY)[__i] += (P[1]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx1+JY)[__i] += (P[1]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx2+JY)[__i] += (P[1]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx6+JY)[__i] += (P[0]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx7+JY)[__i] += (P[0]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx8+JY)[__i] += (P[0]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx9+JY)[__i] += (P[1]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx10+JY)[__i] += (P[1]) * (gz7+JZ)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx11+JY)[__i] += (P[1]) * (gz7+JZ)[__i];

			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx0+JZ)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx1+JZ)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx2+JZ)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx3+JZ)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx4+JZ)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx5+JZ)[__i] += (P[0]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx9+JZ)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx10+JZ)[__i] += (P[1]) * (gz4+JY)[__i];
			for (int __i = 0; __i < (M-JY-JZ); ++__i) (gx11+JZ)[__i] += (P[1]) * (gz4+JY)[__i];

			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];

		} else {
			Real *gs=G+M*6*s_to;
			Real *gs_1=G+M*6*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M;
			Real *gx0=gs,   *gx1=gs+M,   *gx2=gs+2*M,   *gx3=gs+3*M,   *gx4=gs+4*M,   *gx5=gs+5*M;
			Real *g=G1;

			std::fill_n(gs, 6*M, 0);
			set_bounds_x(gz0,gz5,0,0);

			for (int __i = 0; __i < (M-JX); ++__i) (gx1+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4+JX)[__i] += (P[1]) * (gz5)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx5+JX)[__i] += (P[0]) * (gz5)[__i];

			set_bounds_y(gz1,gz4,0,0);

			for (int __i = 0; __i < (M-JY); ++__i) (gx0+JY)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2+JY)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3+JY)[__i] += (P[1]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx4+JY)[__i] += (P[0]) * (gz4)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5+JY)[__i] += (P[1]) * (gz4)[__i];

			set_bounds_z(gz2,gz3,0,0);

			for (int __i = 0; __i < (M-JZ); ++__i) (gx0+JZ)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx1+JZ)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx3+JZ)[__i] += (P[0]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx4+JZ)[__i] += (P[1]) * (gz3)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5+JZ)[__i] += (P[1]) * (gz3)[__i];

			for (int __i = 0; __i < (M-JZ); ++__i) (gx0)[__i] += (P[1]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx1)[__i] += (P[1]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx2)[__i] += (P[0]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx4)[__i] += (P[1]) * (gz2+JZ)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (gx5)[__i] += (P[1]) * (gz2+JZ)[__i];

			for (int __i = 0; __i < (M-JY); ++__i) (gx0)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx1)[__i] += (P[0]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx2)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx3)[__i] += (P[1]) * (gz1+JY)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (gx5)[__i] += (P[1]) * (gz1+JY)[__i];

			for (int __i = 0; __i < (M-JX); ++__i) (gx0)[__i] += (P[0]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx1)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx2)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx3)[__i] += (P[1]) * (gz0+JX)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (gx4)[__i] += (P[1]) * (gz0+JX)[__i];

			for (int k=0; k<6; k++) for (int __i = 0; __i < (M); ++__i) (gs+k*M)[__i] = (gs+k*M)[__i] * (g)[__i];
		}
	} else {
		if (lattice_type==simple_cubic) {
			cout<<"Markov 2, geometry 3, simple_cubic , stencil_full, not implemented" << endl;
		} else {
			cout<<"Markov 2, geometry 3, hexagonal, stencil_full, not implemented" << endl;
		}
	}
}

void LGrad3::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
NAMICS_DBG_THIS(" propagate in LGrad3 " << endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int JX_=JX, JY_=JY;
	int k=sub_box_on;

	std::fill_n(gs, M, 0);
	set_bounds(gs_1);
	if (k>0) {
		JX_=jx[k];
		JY_=jy[k];
	}

	if (!stencil_full) {
		for (int __i = 0; __i < (M-JX_); ++__i) (gs+JX_)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JX_); ++__i) (gs)[__i] += (gs_1+JX_)[__i];
		for (int __i = 0; __i < (M-JY_); ++__i) (gs+JY_)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JY_); ++__i) (gs)[__i] += (gs_1+JY_)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (gs_1+1)[__i];
		if (lattice_type == simple_cubic) {
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (4.0);
		} else {
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (2.0);
		}
		for (int __i = 0; __i < (M-JX_-JY_); ++__i) (gs+JX_+JY_)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_); ++__i) (gs)[__i] += (gs_1+JX_+JY_)[__i];
		for (int __i = 0; __i < (M-JY_-JX_); ++__i) (gs+JY_)[__i] += (gs_1+JX)[__i];
		for (int __i = 0; __i < (M-JY_-JX_); ++__i) (gs+JX)[__i] += (gs_1+JY_)[__i];
		for (int __i = 0; __i < (M-JX_-1); ++__i) (gs+JX_+1)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JX_-1); ++__i) (gs)[__i] += (gs_1+JX_+1)[__i];
		for (int __i = 0; __i < (M-JX_); ++__i) (gs+JX_)[__i] += (gs_1+1)[__i];
		for (int __i = 0; __i < (M-JX_); ++__i) (gs+1)[__i] += (gs_1+JX_)[__i];
		for (int __i = 0; __i < (M-JY_-1); ++__i) (gs+JY_+1)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JX_-1); ++__i) (gs)[__i] += (gs_1+JY_+1)[__i];
		for (int __i = 0; __i < (M-JY_); ++__i) (gs+JY_)[__i] += (gs_1+1)[__i];
		for (int __i = 0; __i < (M-JY_); ++__i) (gs+1)[__i] += (gs_1+JY_)[__i];
		if (lattice_type == simple_cubic) {
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (4.0);
		} else {
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (2.0);
		}
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_+JY_+1)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs)[__i] += (gs_1+JX_+JY_+1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_+JY_)[__i] += (gs_1+1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+1)[__i] += (gs_1+JX_+JY_)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_+1)[__i] += (gs_1+JY_)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JY_)[__i] += (gs_1+JX_+1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JY_+1)[__i] += (gs_1+JX_)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_)[__i] += (gs_1+JY_+1)[__i];
		if (lattice_type == simple_cubic) {
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/152.0);
		} else {
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/56.0);
		}
		for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
	} else {
		if (lattice_type==simple_cubic) {
			for (int __i = 0; __i < (M-JX_); ++__i) (gs+JX_)[__i] += (gs_1)[__i];
			for (int __i = 0; __i < (M-JX_); ++__i) (gs)[__i] += (gs_1+JX_)[__i];
			for (int __i = 0; __i < (M-JY_); ++__i) (gs+JY_)[__i] += (gs_1)[__i];
			for (int __i = 0; __i < (M-JY_); ++__i) (gs)[__i] += (gs_1+JY_)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (gs_1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (gs_1+1)[__i];
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/6.0);
			for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];
		} else { //hexagonal
			if (fjc==1) {
				Real Two=2.0;
				Real C=1.0/40.0;
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] += (Two) * (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (gs_1+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (Two);

				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JX+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs)[__i] += (gs_1+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JX)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY); ++__i) (gs+JY)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX+1)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (gs+1)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (gs+JY+1)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (gs+1)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (gs+JX)[__i] += (gs_1+1)[__i];
				for (int __i = 0; __i < (M-JX-1); ++__i) (gs)[__i] += (gs_1+JX+1)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (gs+JY)[__i] += (gs_1+1)[__i];
				for (int __i = 0; __i < (M-JY-1); ++__i) (gs)[__i] += (gs_1+JY+1)[__i];
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (Two);

				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX+JY+1)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs)[__i] += (gs_1+JX+JY+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX+JY)[__i] += (gs_1+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+1)[__i] += (gs_1+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX+1)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JY)[__i] += (gs_1+JX+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JY+1)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX)[__i] += (gs_1+JY+1)[__i];

				for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (C);
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];


			} else { //hexagonal and fjc=2 ; deze code moet ook werken voor FJC_choices > 5

				for (int block=0; block<4; block++){
					int bk;
					int a,b;
					for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) for (int z=-fjc; z<fjc+1; z++){
						bk=a=b=0;
						if (x==-fjc || x==fjc) bk++;
						if (y==-fjc || y==fjc) bk++;
						if (z==-fjc || z==fjc) bk++;
						if (bk==block) {
							if (x<0) a =-x*JX; else b=x*JX;
							if (y<0) a -=y*JY; else b+=y*JY;
							if (z<0) a -=z*JZ; else b+=z*JZ;
							for (int __i = 0; __i < (M-a-b); ++__i) (gs+a)[__i] += (gs_1+b)[__i];
						}
					}
					if (block !=3) for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (2.0); else for (int __i = 0; __i < (M); ++__i) (gs)[__i] *= (1.0/8.0/((FJC-2)*(FJC-2)*(FJC-2)+3*(FJC-2)*(FJC-2)+3*(FJC-2)+1));

				}
				for (int __i = 0; __i < (M); ++__i) (gs)[__i] = (gs)[__i] * (G1)[__i];


			}
		}
	}
}


bool LGrad3::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
	(void)var_pos;
	(void)n_pos;
	(void)H_p;
NAMICS_DBG_THIS("ReadRange in LGrad3 " << endl);	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	In->split(range,';',set);
	coor.clear();
	block=true; In->split(set[0],',',coor);

	if (coor.size()!=3) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of three: 'x,y,z'" << endl; success=false;}
	else {
		r[0]=In->Get_int(coor[0],0);
		r[1]=In->Get_int(coor[1],0);
		r[2]=In->Get_int(coor[2],0);
	}
	coor.clear(); In->split(set[1],',',coor);

	if (coor.size()!=3) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of three: 'x,y,z'" << endl; success=false;}
	else {
		r[3]=In->Get_int(coor[0],0);
		r[4]=In->Get_int(coor[1],0);
		r[5]=In->Get_int(coor[2],0);
	}
	if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	if (r[2] > r[5]) {cout << "In mon " + seg_name+ ", for 'pos 1', the z-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	return success;
}

bool LGrad3::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
NAMICS_DBG_THIS("ReadRangeFile in LGrad3 " << endl);	if (fjc>1) {
		cout << "Rangefile is not implemented for FJC-choices >3; contact FL. " << endl;
		return false;
	}

	bool success=true;
	string content;
	vector<string> lines;
	vector<string> sub;
	vector<string> xyz;
	string Infilename=In->name;
	In->split(Infilename,'.',sub);

	int length;
	int length_xyz;
	int px,py,pz,p_i,x,y,z;
	int i=0;
	if (!range_reader->ReadSanitizedFile(sub[0].append(".").append(filename),content)) {
		success=false;
		return success;
	}

	In->split(content,'#',lines);
	length = lines.size();
	if (length == MX*MY*MZ) { //expect to read 'mask file';
		i=0;
		if (n_pos==0) {
			while (i<length){
				if (In->Get_int(lines[i],0)==1) n_pos++;
				i++;
			};
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
		} else {
			i=0; p_i=0;
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++) {
				if (In->Get_int(lines[i],0)==1) {H_p[p_i]=x*JX+y*JY+fjc-1+z; p_i++;}
				i++;
			}
		}
	} else { //expect to read x,y,z
		px=0,py=0,pz=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=3) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'triple coordinate' structure 'x,y,z' was not found. " << endl;  success = false;
				} else {
					px=In->Get_int(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
					py=In->Get_int(xyz[1],0);
					if (py < 1 || py > MY) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MY << endl; success =false;}
					pz=In->Get_int(xyz[2],0);
					if (pz < 1 || pz > MZ) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MZ << endl; success =false;}
				}
				H_p[i]=px*JX+py*JY+fjc-1+pz;
				i++;
			}
		}
	}
	return success;
}

bool LGrad3::FillMask(Real* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
	bool success=true;
	bool readfile=false;
	int length=0;
	int length_px = px.size();

	vector<string> lines;
	int p;
	if (px.size()==0) {
		readfile=true;
		string content;
		success=range_reader->ReadSanitizedFile(filename,content);
		if (success) {
			In->split(content,'#',lines);
			length = lines.size();
		}
	}

	if (readfile) {
		if (MX*MY*MZ!=length) {
			success=false; cout <<"inputfile for filling delta_range has not expected length in x,y,z-directions" << endl;
		} else {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++)
			for (int z=1; z<MZ+1; z++) Mask[x*JX + y*JY + z]=In->Get_int(lines[x*JX + y*JY + z],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			p=py[i]; if (p<1 || p>MY) {success=false; cout<<" y-value in delta_range out of bounds; " << endl; }
			p=pz[i]; if (p<1 || p>MZ) {success=false; cout<<" z-value in delta_range out of bounds; " << endl; }
			if (success) Mask[px[i]*JX + py[i]*JY + fjc-1+ pz[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad3::CreateMASK(Real* H_MASK, int* r, int* H_P, int n_pos, bool block) {
NAMICS_DBG_THIS("CreateMask for LGrad3 " + name << endl);	bool success=true;
	std::fill_n(H_MASK, M, static_cast<Real>(0));
	// Build mask from either a block in r=[x1,y1,z1,x2,y2,z2] or list of indices in H_P.
	if (block) {
		// mark all (x,y,z) in [x1, x2] x [y1, y2] x [z1, z2].
		for (int x=r[0]; x<r[3]+1; x++)
		for (int y=r[1]; y<r[4]+1; y++)
		for (int z=r[2]; z<r[5]+1; z++)
			H_MASK[P(x,y,z)]=1;
	} else {
		// mark n_pos linear indices from H_P.
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad3::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	(result) = 0; for (int __i = 0; __i < (M); ++__i) (result) += (phi)[__i] * (L)[__i];
	return result;
}

void LGrad3::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);

	std::fill_n(EE, M, 0);
	for (int __i = 0; __i < (M-2); ++__i) (EE+1)[__i] += pow((psi)[__i]-(psi+1)[__i],2) + pow((psi+1)[__i]-(psi+2)[__i],2);
	for (int __i = 0; __i < (M-2*JX); ++__i) (EE+JX)[__i] += pow((psi)[__i]-(psi+JX)[__i],2) + pow((psi+JX)[__i]-(psi+2*JX)[__i],2);
	for (int __i = 0; __i < (M-2*JY); ++__i) (EE+JY)[__i] += pow((psi)[__i]-(psi+JY)[__i],2) + pow((psi+JY)[__i]-(psi+2*JY)[__i],2);
	for (int __i = 0; __i < (M); ++__i) (EE)[__i] *= (pf);


}


void LGrad3::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, Real* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	(void)grad_epsilon;
	int x, y, z;

	Real epsZplus, epsZmin, epsXplus, epsXmin, epsYplus, epsYmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {


	for (x=1; x<MX+1; x++) {
		for (y=1; y<MY+1; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ+1; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
				}
			}
		}
	}
	for (x=MX; x>0; x--) {
		for (y=MY; y>0; y--) {
			epsZmin=eps[x*JX+y*JY+MZ+1]+eps[x*JX+y*JY+MZ];
			for (z=MZ; z>0; z--) {
				epsZplus=epsZmin;
				epsZmin=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z-1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
					g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
				}
			}
		}
	}





   } else { //fixedPsi0 is true

	for (x=1; x<MX+1; x++) {
		for (y=1; y<MY+1; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ+1; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0)
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					    epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
				           epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					    C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
			}
		}
	}
	for (x=MX; x>0; x--) {
		for (y=MY; y>0; y--) {
			epsZmin=eps[x*JX+y*JY+MZ+1]+eps[x*JX+y*JY+MZ];
			for (z=MZ; z>0; z--) {
				epsZplus=epsZmin;
				epsZmin=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z-1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
					g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
				}
			}
		}
	}




   }
}


void LGrad3::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, Real* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	(void)grad_epsilon;
	int z, x, y;
	Real epsXplus,epsXmin,epsYplus,epsYmin,epsZplus,epsZmin;

	Real C = -e*e/(eps0*k_BT*bond_length);

	for (x=1; x<MX; x++) {
		for (y=1; y<MY; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==1) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					                   epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
							     epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]-
							     (epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus)*psi[x*JX+y*JY+z])/C;
					g[x*JX+y*JY+z]=-q[x*JX+y*JY+z];
				}
			}
		}
	}





}

void LGrad3::set_bounds_x(Real* X,Real* Y,int shifty,int shiftz){
NAMICS_DBG_THIS("set_bounds_x (shift in y,z ) in LGrad3 " << endl);	int y,z;
	int k=0;
	if (BX1>BXM) {
		set_bounds_x(X,shifty,shiftz);
		set_bounds_x(Y,shifty,shiftz);
	} else {

		if (fjc==1) {
			 for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
				X[0+        y*JY+z*JZ] = Y[BX1*JX+(y+shifty)*JY+(z+shiftz)*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = Y[BXM*JX+(y-shifty)*JY+(z-shiftz)*JZ];
				Y[0        +y*JY+z*JZ] = X[BX1*JX+(y+shifty)*JY+(z+shiftz)*JZ];
				Y[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+(y-shifty)*JY+(z-shiftz)*JZ];
			}
		} else {
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) {
					X[k*JX+y*JY+z*JZ] = Y[B_X1[k]*JX+(y+shifty)*JY+(z+shiftz)*JZ];
					Y[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+(y+shifty)*JY+(z+shiftz)*JZ];

				}
				for (k=0; k<fjc; k++) {
					X[(MX+fjc+k)*JX+y*JY+z*JZ] = Y[B_XM[k]*JX+(y-shifty)*JY+(z-shiftz)*JZ];
					Y[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+(y-shifty)*JY+(z-shiftz)*JZ];
				}
			}
		}
	}
}
void LGrad3::set_bounds_y(Real* X,Real* Y,int shiftx, int shiftz){
NAMICS_DBG_THIS("set_bounds_y (shift in x,z) in LGrad3 " << endl);	int x,z;
	int k=0;
	if (BY1>BYM) {
		set_bounds_y(X,shiftx,shiftz);
		set_bounds_y(Y,shiftx,shiftz);
	} else {

		if (fjc==1) {
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0+        z*JZ] = Y[(x+shiftx)*JX+BY1*JY+(z+shiftz)*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = Y[(x-shiftx)*JX+BYM*JY+(z-shiftz)*JZ];
				Y[x*JX+0+        z*JZ] = X[(x+shiftx)*JX+BY1*JY+(z+shiftz)*JZ];
				Y[x*JX+(MY+1)*JY+z*JZ] = X[(x-shiftx)*JX+BYM*JY+(z-shiftz)*JZ];
			}
		} else {
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) {
					X[x*JX+k*JY+z*JZ] = Y[(x+shiftx)*JX+B_Y1[k]*JY+(z+shiftz)*JZ];
					Y[x*JX+k*JY+z*JZ] = X[(x+shiftx)*JX+B_Y1[k]*JY+(z+shiftz)*JZ];
				}
				for (k=0; k<fjc; k++) {
					X[x*JX+(MY+fjc+k)*JY+z*JZ] = Y[(x-shiftx)*JX+B_YM[k]*JY+(z-shiftz)*JZ];
					Y[x*JX+(MY+fjc+k)*JY+z*JZ] = X[(x-shiftx)*JX+B_YM[k]*JY+(z-shiftz)*JZ];
				}
			}
		}
	}

}

void LGrad3::set_bounds_z(Real* X,Real* Y,int shiftx,int shifty){
NAMICS_DBG_THIS("set_bounds_z shift (x,y) in LGrad3 " << endl);	int x,y;
	int k=0;
	if (BZ1>BZM) { //periodic
		set_bounds_z(X,shiftx,shifty);
		set_bounds_z(Y,shiftx,shifty);
	} else {
		if (fjc==1) {
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) {
				X[x*JX+y*JY+   0] = Y[(x+shiftx)*JX+(y+shifty)*JY+BZ1];
				X[x*JX+y*JY+MZ+1] = Y[(x-shiftx)*JX+(y-shifty)*JY+BZM];
				Y[x*JX+y*JY+   0] = X[(x+shiftx)*JX+(y+shifty)*JY+BZ1];
				Y[x*JX+y*JY+MZ+1] = X[(x-shiftx)*JX+(y-shifty)*JY+BZM];
			}
		} else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) {
				for (k=0; k<fjc; k++) {
					X[x*JX+y*JY+k] = Y[(x+shiftx)*JX+(y+shifty)*JY+B_Z1[k]];
					Y[x*JX+y*JY+k] = X[(x+shiftx)*JX+(y+shifty)*JY+B_Z1[k]];
				}
				for (k=0; k<fjc; k++) {
					X[x*JX+y*JY+MZ+fjc+k]  = Y[(x-shiftx)*JX+(y-shifty)*JY+B_ZM[k]];
					Y[x*JX+y*JY+MZ+fjc+k]  = X[(x-shiftx)*JX+(y-shifty)*JY+B_ZM[k]];
				}
			}
		}
	}
}

void LGrad3::set_bounds_x(Real* X, int shifty, int shiftz){
NAMICS_DBG_THIS("set_bounds_x (shift yz) in LGrad3 " << endl);	int y,z;
	int k=0;
	if (fjc==1) {
		 for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
			X[0        +y*JY+z*JZ] = X[BX1*JX+(y+shifty)*JY+(z+shiftz)*JZ];
			X[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+(y-shifty)*JY+(z-shiftz)*JZ];
		}
	} else {//shifts for hexagonal lattice not put in the fjc>1 case...
		for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
			for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+(y+shifty)*JY+(z+shiftz)*JZ];
			for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+(y-shifty)*JY+(z-shiftz)*JZ];
		}
	}
}

void LGrad3::set_bounds_y(Real* X,int shiftx, int shiftz){
NAMICS_DBG_THIS("set_bounds_y (shift x,z) in LGrad3 " << endl);	int x,z;
	int k=0;
	if (fjc==1) {
		for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
			X[x*JX+    0    +z*JZ] = X[(x+shiftx)*JX+BY1*JY+(z+shiftz)*JZ];
			X[x*JX+(MY+1)*JY+z*JZ] = X[(x-shiftx)*JX+BYM*JY+(z-shiftz)*JZ];
		}
	} else {
		for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
			for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[(x+shiftx)*JX+B_Y1[k]*JY+(z+shiftz)*JZ];
			for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[(x-shiftx)*JX+B_YM[k]*JY+(z-shiftz)*JZ];
		}
	}
}

void LGrad3::set_bounds_z(Real* X,int shiftx, int shifty){
NAMICS_DBG_THIS("set_bounds_z (shift xy) in LGrad3 " << endl);	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) {
			X[x*JX+y*JY+   0] = X[(x+shiftx)*JX+(y+shifty)*JY+BZ1];
			X[x*JX+y*JY+MZ+1] = X[(x-shiftx)*JX+(y-shifty)*JY+BZM];
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
			for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[(x+shiftx)*JX+(y+shifty)*JY+B_Z1[k]];
			for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[(x-shiftx)*JX+(y-shifty)*JY+B_ZM[k]];
		}
	}
}

void LGrad3::remove_bounds(Real *X){
NAMICS_DBG_THIS("remove_bounds in LGrad3 " << endl);	int x,y,z;
	int k;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			RemoveBoundaries(std::span<Real>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) RemoveBoundaries(std::span<Real>(X, static_cast<size_t>(M)),JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=0; x<MX+2*fjc; x++) for (y=0; y<MY+2*fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = 0;
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = 0;
			}
			for (y=0; y<MY+2*fjc; y++) for (z=0; z<MZ+2*fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = 0;
			}
			for (z=0; z<MZ+2*fjc; z++) for (x=0; x<MX+2*fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = 0;
			}
		}
	}
}

void LGrad3::set_bounds(Real* X){
NAMICS_DBG_THIS("set_bounds in LGrad3 " << endl);	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(std::span<Real>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++){
				X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
				X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
			}
			for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
				X[0        +y*JY+z*JZ] = X[BX1*JX+y*JY+z*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+y*JY+z*JZ];
			}
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
			}

			x=0; {
				for (y=1; y<MY+1; y++){
					X[y*JY+0]     = X[y*JY+BZ1];
					X[y*JY+MZ+1]  = X[y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[0        +z*JZ] = X[BY1*JY+z*JZ];
					X[(MY+1)*JY+z*JZ] = X[BYM*JY+z*JZ];
				}
			}
			x=MX+1; {
				for (y=1; y<MY+1; y++){
					X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
					X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
					X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
				}
			}
			y=0; {
				for (x=1; x<MX+1; x++){
					X[x*JX        +0] = X[x*JX+BZ1*JZ];
					X[x*JX+(MZ+1)*JZ] = X[x*JX+BZM*JZ];
				}
			}
			y=MY+1; {
				for (x=1; x<MX+1; x++){
					X[x*JX+y*JY        +0] = X[x*JX+y*JY+BZ1*JZ];
					X[x*JX+y*JY+(MZ+1)*JZ] = X[x*JX+y*JY+BZM*JZ];
				}
			}

			X[0        +0        +0]        =X[BX1*JX+BY1*JY+BZ1*JZ];
			X[(MX+1)*JX+0        +0]        =X[BXM*JX+BY1*JY+BZ1*JZ];
			X[0        +(MY+1)*JY+0]        =X[BX1*JX+BYM*JY+BZ1*JZ];
			X[0        +0        +(MZ+1)*JZ]=X[BX1*JX+BY1*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+0]        =X[BXM*JX+BYM*JY+BZ1*JZ];
			X[(MX+1)*JX+ 0       +(MZ+1)*JZ]=X[BXM*JX+BY1*JY+BZM*JZ];
			X[0        +(MY+1)*JY+(MZ+1)*JZ]=X[BX1*JX+BYM*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+(MZ+1)*JZ]=X[BXM*JX+BYM*JY+BZM*JZ];

		} else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[x*JX+y*JY+B_Z1[k]];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[x*JX+B_Y1[k]*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+B_YM[k]*JY+z*JZ];
			}

			for (x=0; x<fjc; x++ ) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) X[x*JX+k*JY         +z]  = X[x*JX+B_Y1[k]*JY+z*JZ];
					for (int k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z*JZ];
				}
			}
			for (x=MX+fjc; x<MX+2*fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) X[x*JX+k*JY        +z]  = X[x*JX+B_Y1[k]*JY+z];
					for (int k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z];
				}
			}
			for (y=0; y<fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) X[x*JX +y*JY       +k]    = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX +y*JY +(MZ+fjc+k)] = X[x*JX+y*JY+B_ZM[k]];
				}
			}
			for (y=MY+fjc; y<MY+2*fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY        +k]     = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+(MZ+fjc+k)*JZ] = X[x*JX+y*JY+B_ZM[k]];
				}
			}

			for (int k=0; k<fjc; k++) for (int l=0; l<fjc; l++) for (int m=0; m<fjc; m++) {
				X[k*JX         +l*JY         +m*JZ]         =X[B_X1[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +m*JZ]         =X[B_XM[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+m*JZ]         =X[B_X1[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +l*JY         +(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+m*JZ]         =X[B_XM[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
			}
		}
	}
}

void LGrad3::set_M_bounds(Real* X){
if (!debug) cout <<"set_bounds (M) in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(std::span<Real>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {
			//SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++){
				X[x*JX+y*JY+0]     = X[x*JX+y*JY+1];
				X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+MZ];
			}
			for (y=1; y<MY+1; y++) for (z=fjc; z<MZ+1; z++)  {
				X[0        +y*JY+z*JZ] = X[1*JX+y*JY+z*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = X[MX*JX+y*JY+z*JZ];
			}
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0        +z*JZ] = X[x*JX+1*JY+z*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+MY*JY+z*JZ];
			}

		} else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[x*JX+y*JY+2*fjc-1-k];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+MZ+fjc-k-1];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[(2*fjc-1-k)*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[(MX+fjc-k-1)*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){ //corners?
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[x*JX+(2*fjc-1-k)*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+(MY+fjc-k-1)*JY+z*JZ];
			}
		}
	}
}

void LGrad3::remove_bounds(int *X){
if (!debug) cout <<"remove_bounds (int) in LGrad3 " << endl;
	int x,y,z;
	int k;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			RemoveBoundaries(std::span<int>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) RemoveBoundaries(std::span<int>(X, static_cast<size_t>(M)),JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=0; x<MX+2*fjc; x++) for (y=0; y<MY+2*fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = 0;
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = 0;
			}
			for (y=0; y<MY+2*fjc; y++) for (z=0; z<MZ+2*fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = 0;
			}
			for (z=0; z<MZ+2*fjc; z++) for (x=0; x<MX+2*fjc; x++){//corners?
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = 0;
			}
		}
	}
}

void LGrad3::set_bounds(int* X){
if (!debug) cout <<"set_bounds (int) in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(std::span<int>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {
			//SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

			for (x=1; x<MX+1; x++) {
				for (y=1; y<MY+1; y++){
					X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
					X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++) {
					X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
					X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
				}
			}
			for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
				X[0        +y*JY+z*JZ] = X[BX1*JX+y*JY+z*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+y*JY+z*JZ];
			}


			for (int x=0; x<MX+2; x+=MX+1) {
				for (y=1; y<MY+1; y++){
					X[y*JY+0]     = X[y*JY+BZ1];
					X[y*JY+MZ+1]  = X[y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[0        +z*JZ] = X[BY1*JY+z*JZ];
					X[(MY+1)*JY+z*JZ] = X[BYM*JY+z*JZ];
				}
			}

			for (y=0; y<MY+2; y+=MY+1) {
				for (int x=1; x<MX+1; x++){
					X[x*JX        +0] = X[x*JX+BZ1*JZ];
					X[x*JX+(MZ+1)*JZ] = X[x*JX+BZM*JZ];
				}
			}

			X[0        +0        +0]        =X[BX1*JX+BY1*JY+BZ1*JZ];
			X[(MX+1)*JX+0        +0]        =X[BXM*JX+BY1*JY+BZ1*JZ];
			X[0        +(MY+1)*JY+0]        =X[BX1*JX+BYM*JY+BZ1*JZ];
			X[0        +0        +(MZ+1)*JZ]=X[BX1*JX+BY1*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+0]        =X[BXM*JX+BYM*JY+BZ1*JZ];
			X[(MX+1)*JX+ 0       +(MZ+1)*JZ]=X[BXM*JX+BY1*JY+BZM*JZ];
			X[0        +(MY+1)*JY+(MZ+1)*JZ]=X[BX1*JX+BYM*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+(MZ+1)*JZ]=X[BXM*JX+BYM*JY+BZM*JZ];


		}else {
			for (x=fjc; x<MX+fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){
					for (k=0; k<fjc; k++) {
						X[x*JX+y*JY+k] = X[x*JX+y*JY+B_Z1[k]];
						X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
					}
				}
				for (z=fjc; z<MZ+fjc; z++) {
					for (k=0; k<fjc; k++) {
						X[x*JX+k*JY+z*JZ] = X[x*JX+B_Y1[k]*JY+z*JZ];
						X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+B_YM[k]*JY+z*JZ];
					}
				}
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) {
					X[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+y*JY+z*JZ];
					X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+y*JY+z*JZ];
				}
			}

			for (x=0; x<fjc; x++ ) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) {
						X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
						X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
					}
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) {
						X[x*JX+k*JY         +z]  = X[x*JX+B_Y1[k]*JY+z*JZ];
						X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z*JZ];
					}
				}
			}
			for (x=MX+fjc; x<MX+2*fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) {
						X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
						X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
					}
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) {
						X[x*JX+k*JY        +z]  = X[x*JX+B_Y1[k]*JY+z];
						X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z];
					}
				}
			}
			for (y=0; y<fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) {
						X[x*JX +y*JY       +k]    = X[x*JX+y*JY+B_Z1[k]];
						X[x*JX +y*JY +(MZ+fjc+k)] = X[x*JX+y*JY+B_ZM[k]];
					}
				}
			}
			for (y=MY+fjc; y<MY+2*fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) {
						X[x*JX+y*JY        +k]     = X[x*JX+y*JY+B_Z1[k]];
						X[x*JX+y*JY+(MZ+fjc+k)*JZ] = X[x*JX+y*JY+B_ZM[k]];
					}
				}
			}

			for (int k=0; k<fjc; k++) for (int l=0; l<fjc; l++) for (int m=0; m<fjc; m++) {
				X[k*JX         +l*JY         +m*JZ]         =X[B_X1[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +m*JZ]         =X[B_XM[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+m*JZ]         =X[B_X1[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +l*JY         +(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+m*JZ]         =X[B_XM[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
			}
		}
	}
}

Real LGrad3::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	int size;

	if (Markov==2) {
		if (lattice_type==simple_cubic) size =6;
		if (lattice_type==hexagonal) size =12;

		for (int k=0; k<size; k++) GN+=WeightedSum(G+k*M);
		GN/=(1.0*size);
	} else GN=WeightedSum(G);
	return GN;

}
void LGrad3::AddPhiS(Real* phi,Real* Gf,Real* Gb,int Markov, int M){
	Real one=1.0;
	if (Markov ==2) {
		if (lattice_type==hexagonal) {
			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (one/12.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
		} else {
			for (int k=0; k<6; k++) for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (one/6.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
		}
	} else for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (Gf)[__i] * (Gb)[__i];

}

void LGrad3::AddPhiS(Real* phi,Real* Gf,Real* Gb,Real degeneracy,int Markov, int M){
	if (Markov ==2) {
		if (lattice_type==hexagonal) {
			for (int k=0; k<12; k++) for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/12.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
		} else {
			for (int k=0; k<6; k++) for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy/6.0) * (Gf+k*M)[__i] * (Gb+k*M)[__i];
		}
	} else for (int __i = 0; __i < (M); ++__i) (phi)[__i] += (degeneracy) * (Gf)[__i] * (Gb)[__i];

}


void LGrad3::AddPhiS(Real* phi,Real* Gf,Real* Gb, Real* G1, Real norm, int Markov, int M){
	(void)M;
	(void)Markov;
	(void)norm;
	(void)G1;
	(void)Gb;
	(void)Gf;
	(void)phi;
	cout << "Composition phi Alias not implemented in LGrad3 " << endl;
}

void LGrad3::Initiate(Real* G,Real* Gz,int Markov, int M){
	int size;
	if (Markov==2){
		if (lattice_type==simple_cubic) size =6;
		if (lattice_type==hexagonal) {
			size =12;
		}
		for (int k=0; k<size; k++) std::copy_n(Gz, M, G+k*M);
	} else std::copy_n(Gz, M, G);
}

void LGrad3::Terminate(Real* Gz,Real* G,int Markov, int M){
	int size;
	if (Markov==2){
		std::fill_n(Gz, M, 0);
		if (lattice_type==simple_cubic) size =6;
		if (lattice_type==hexagonal) {
			size =12;
		}
		for (int k=0; k<size; k++) for (int __i = 0; __i < (M); ++__i) (Gz)[__i] += (G+k*M)[__i];
		for (int __i = 0; __i < (M); ++__i) (Gz)[__i] *= (1.0/size);
	} else std::copy_n(G, M, Gz);
}

bool LGrad3:: PutMask(Real* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	bool success=true;
	int length =px.size();
	int X,Y,Z;
	std::fill_n(MASK, M, 0);

	for (int i =0; i<length; i++) {
		int xx,yy,zz;
		xx=px[i]; yy=py[i]; zz=pz[i];
		for (int x=xx-R; x<xx+R+1; x++)
		for (int y=yy-R; y<yy+R+1; y++)
		for (int z=zz-R; z<zz+R+1; z++) {
			if ((xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z) <=R*R) {
				X=x; Y=y; Z=z;
				if (x<1) {if (BX1==1) X=0; else X+=MX;}
				if (y<1) {if (BY1==1) Y=0; else Y+=MY;}
				if (z<1) {if (BZ1==1) Z=0; else Z+=MZ;}
				if (x>MX) {if (BXM==MX) X=MX; else X-=MX;}
				if (y>MY) {if (BYM==MY) Y=MY; else Y-=MY;}
				if (z>MZ) {if (BZM==MZ) Z=MZ; else Z-=MZ;}
				MASK[P(X,Y,Z)]=1;
			}
		}
		//}
	}

	return success;
}


Real LGrad3::DphiDt(Real* g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB, Real B_A, Real B_B) {
	(void)B_B;
	(void)B_A;
	(void)alphaB;
	(void)alphaA;
	(void)phiB;
	(void)phiA;
	(void)B_phitot;
	(void)g;
	cout <<"Grad3 DphiDt not implemented yet " << endl;
	return 0;
}

Real LGrad3::MomentPlanar(Real* X,int n,Real Z0){
	(void)Z0;
	(void)n;
	(void)X;
	cout <<"MomentPlanar not implemented; kJ0 or kbar may be wrong. " << endl;
	return 0;
}
