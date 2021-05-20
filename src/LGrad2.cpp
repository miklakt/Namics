#include <iostream> 
#include <string> 
#include "lattice.h"
#include "LGrad2.h" 

LGrad2::LGrad2(vector<Input*> In_,string name_): Lattice(In_,name_) {}

LGrad2::~LGrad2() {
if (debug) cout <<"LGrad2 destructor " << endl;
}


void LGrad2:: ComputeLambdas() {
	Real r, VL, LS;
	Real rlow, rhigh;


	if (fcc_sites){
			if (geometry=="cylindrical"){
				for (int x=1; x<MX+1; x++)
				for (int y=1; y<MY+1; y++) {
					r=offset_first_layer + 1.0*x;
					L[P(x,y)]=PIE*(pow(r,2)-pow(r-1,2));
					fcc_lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]/3.0;
					fcc_lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]/3.0;
					fcc_lambda0[P(x,y)]=1.0-2.0/3.0;
				}
			}
	}
	if (fjc==1) {
		if (geometry=="planar") {
			for (int i=0; i<M; i++) L[i]=1;
		}
		if (geometry=="cylindrical") {
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
		}
	}
	if (fjc==2) {
		if (geometry=="planar") {
			for (int i=0; i<M; i++) L[i]=1.0/fjc;
		}
		if (geometry == "cylindrical") {
			for (int y=1-fjc; y<MY+fjc; y++) {
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
	}
	for (int k=0; k<FJC; k++) for (int y=1-fjc; y<MY+fjc; y++) for (int x=MX+1; x< MX+fjc; x++)
		LAMBDA[P(x,y)+k*M]=LAMBDA[P(2*MX-x+1,y)+(FJC-k-1)*M];
}

bool LGrad2::PutM() {
if (debug) cout << "PutM in LGrad2 " << endl;
	bool success=true;

	if (geometry=="cylindrical")
		volume = MY*PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));
	else volume = MX*MY;
	JX=MY+2*fjc; JY=1; JZ=0; M=(MX+2*fjc)*(MY+2*fjc);

	Accesible_volume=volume;
	return success;
}

void LGrad2::TimesL(Real* X){
if (debug) cout << "TimesL in LGrad2 " << endl;
	if (geometry!="planar") Times(X,X,L,M);
}

void LGrad2::DivL(Real* X){
if (debug) cout << "DivL in LGrad2 " << endl;
	if (geometry!="planar") Div(X,L,M);
}

Real LGrad2:: Moment(Real* X,Real Xb, int n) {
if (debug) cout << "Moment in LGrad2 " << endl;
	Real Result=0;
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
		cout <<"Moment analysis not yet implemented in LGrad" << endl;
	}
	return Result/fjc;
}

Real LGrad2::WeightedSum(Real* X){
if (debug) cout << "weighted sum in LGrad2 " << endl;
	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		Sum(sum,X,M); sum = sum/(fjc*fjc);
	} else	{
		Dot(sum,X,L,M);
	}
	return sum;
}

void LGrad2::vtk(string filename, Real* X, string id,bool writebounds) {
if (debug) cout << "vtk in LGrad2 " << endl;
	FILE *fp;
	int i;
	fp = fopen(filename.c_str(),"w+");
	fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MY,MX,1);
	if (writebounds) {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",(MX+2*fjc)*(MY+2*fjc));
	} else {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY);
	}
	fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
	if (writebounds) for (i=0; i<M; i++) fprintf(fp,"%f\n",X[i]);
	for (int x=1; x<MX+1; x++)
	for (int y=1; y<MY+1; y++)
	fprintf(fp,"%f\n",X[P(x,y)]);
	fclose(fp);
}

void LGrad2::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds){
if (debug) cout <<"PutProfiles in LGrad2 " << endl;
	int x,y,i;
	int length=X.size();
	int a;
	if (writebounds) a=fjc; else a = 0;
	for (x=1-a; x<MX+1+a; x++)
	for (y=1-a; y<MY+1+a; y++){
		fprintf(pf,"%e\t%e\t",offset_first_layer+1.0*x/fjc-1/(2.0*fjc),1.0*y/fjc-1/(2.0*fjc));
		for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][P(x,y)]);
		fprintf(pf,"\n");
	}
}


void LGrad2::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LGrad2 " << endl;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);set_bounds(X);

	if (fcc_sites) {
			if (geometry=="planar") {
				YplusisCtimesX(X_side,X,     1.0/9.0,M);
				YplusisCtimesX(X_side+1,X,   1.0/9.0,M-1);
				YplusisCtimesX(X_side,X+1,   1.0/9.0,M-1);
				YplusisCtimesX(X_side+JX,X,  1.0/9.0,M-JX);
				YplusisCtimesX(X_side,X+JX,  1.0/9.0,M-JX);
				YplusisCtimesX(X_side+JX+1,X,1.0/9.0,M-JX-1);
				YplusisCtimesX(X_side+JX,X+1,1.0/9.0,M-JX);
				YplusisCtimesX(X_side+1,X+JX,1.0/9.0,M-JX);
				YplusisCtimesX(X_side,X+JX+1,1.0/9.0,M-JX-1);
			} else {
				YplusisCtimesX(X_side,X,1.0/3.0,M);
				AddTimes(X_side+JX,X,fcc_lambda_1+JX,M-JX);
				AddTimes(X_side,X+JX,fcc_lambda1,M-JX);
				YplusisCtimesX(X_side+1,X,1.0/3.0,M-1);
				YplusisCtimesX(X_side,X+1,1.0/3.0,M-1);
				AddTimes(X_side+JX+1,X,fcc_lambda_1+JX+1,M-JX-1);
				AddTimes(X_side+JX,X+1,fcc_lambda_1+JX,M-JX);
				AddTimes(X_side+1,X+JX,fcc_lambda1+1,M-JX);
				AddTimes(X_side,X+JX+1,fcc_lambda1,M-JX-1);
				Norm(X_side,1.0/3.0,M);
			}
	} else {
		if (fjc==1) {
			if (geometry=="planar") {
				if (lattice_type=="simple_cubic") {
					YplusisCtimesX(X_side,X,    16.0/36.0,M);
					YplusisCtimesX(X_side+1,X,   4.0/36.0,M-1);
					YplusisCtimesX(X_side,X+1,   4.0/36.0,M-1);
					YplusisCtimesX(X_side+JX,X,  4.0/36.0,M-JX);
					YplusisCtimesX(X_side,X+JX,  4.0/36.0,M-JX);
					YplusisCtimesX(X_side+JX+1,X,1.0/36.0,M-JX-1);
					YplusisCtimesX(X_side+JX,X+1,1.0/36.0,M-JX);
					YplusisCtimesX(X_side+1,X+JX,1.0/36.0,M-JX);
					YplusisCtimesX(X_side,X+JX+1,1.0/36.0,M-JX-1);
							//6point stencil
							//Add(X_side+JX,X,M-JX);
							//Add(X_side,X+JX,M-JX);
							//Add(X_side+1,X,M-1);
							//Add(X_side,X+1,M-1);
							//Norm(X_side,1.0/2.0,M);
							//Add(X_side,X,M);
							//Norm(X_side,1.0/3.0,M);
				} else {
					YplusisCtimesX(X_side,X,    12.0/48.0,M);
					YplusisCtimesX(X_side+1,X,   6.0/48.0,M-1);
					YplusisCtimesX(X_side,X+1,   6.0/48.0,M-1);
					YplusisCtimesX(X_side+JX,X,  6.0/48.0,M-JX);
					YplusisCtimesX(X_side,X+JX,  6.0/48.0,M-JX);
					YplusisCtimesX(X_side+JX+1,X,3.0/48.0,M-JX-1);
					YplusisCtimesX(X_side+JX,X+1,3.0/48.0,M-JX);
					YplusisCtimesX(X_side+1,X+JX,3.0/48.0,M-JX);
					YplusisCtimesX(X_side,X+JX+1,3.0/48.0,M-JX-1);
				}
			} else {
				if (lattice_type =="simple_cubic") {
					YplusisCtimesX(X_side,X,4.0/6.0,M);
					AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
					AddTimes(X_side,X+JX,lambda1,M-JX);
					YplusisCtimesX(X_side+1,X,1.0/6.0,M-1);
					YplusisCtimesX(X_side,X+1,1.0/6.0,M-1);
					Norm(X_side,4.0,M);
					AddTimes(X_side+JX+1,X,lambda_1+JX+1,M-JX-1);
					AddTimes(X_side+JX,X+1,lambda_1+JX,M-JX);
					AddTimes(X_side+1,X+JX,lambda1+1,M-JX);
					AddTimes(X_side,X+JX+1,lambda1,M-JX-1);
					Norm(X_side,1.0/6.0,M);
				} else {
					YplusisCtimesX(X_side,X,2.0/4.0,M);
					AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
					AddTimes(X_side,X+JX,lambda1,M-JX);
					YplusisCtimesX(X_side+1,X,1.0/4.0,M-1);
					YplusisCtimesX(X_side,X+1,1.0/4.0,M-1);
					Norm(X_side,2.0,M);
					AddTimes(X_side+JX+1,X,lambda_1+JX+1,M-JX-1);
					AddTimes(X_side+JX,X+1,lambda_1+JX,M-JX);
					AddTimes(X_side+1,X+JX,lambda1+1,M-JX);
					AddTimes(X_side,X+JX+1,lambda1,M-JX-1);
					Norm(X_side,3.0/12.0,M);
				}
			}
		}
		if (fjc==2) {
			if (geometry=="planar") {
				Add(X_side,X,M);
				Add(X_side+JX,X,M-JX);
				Add(X_side,X+JX,M-JX);
				Add(X_side+1,X,M-1);
				Add(X_side,X+1,M-1);
				Add(X_side+JX+1,X,M-JX-1);
				Add(X_side,X+JX+1,M-JX-1);
				Add(X_side+1,X+JX,M-JX);
				Add(X_side+JX,X+1,M-JX);
				Norm(X_side,2.0,M);
				Add(X_side+2*JX,X,M-2*JX);
				Add(X_side,X+2*JX,M-2*JX);
				Add(X_side+2*JX,X+1,M-2*JX);
				Add(X_side+2*JX+1,X,M-2*JX-1);
				Add(X_side+1,X+2*JX,M-2*JX);
				Add(X_side,X+2*JX+1,M-2*JX-1);
				Add(X_side+JX+2,X,M-JX-2);
				Add(X_side+JX,X+2,M-JX);
				Add(X_side,X+JX+2,M-JX-2);
				Add(X_side+2,X+JX,M-JX);
				Add(X_side+2,X,M-2);
				Add(X_side,X+2,M-2);
				Norm(X_side,2.0,M);
				Add(X_side+2*JX+2,X,M-2*JX-2);
				Add(X_side,X+2*JX+2,M-2*JX-2);
				Add(X_side+2*JX,X+2,M-2*JX);
				Add(X_side+2,X+2*JX,M-2*JX);
				Norm(X_side,1.0/64.0,M);
			} else {
				Add(X_side,X,M);
				Add(X_side+1,X,M-1);
				Add(X_side,X+1,M-1);
				Norm(X_side,1.0/4.0,M);
				AddTimes(X_side+JX,X,LAMBDA+M+JX,M-JX);
				AddTimes(X_side+JX+1,X,LAMBDA+M+JX+1,M-JX-1);
				AddTimes(X_side+JX,X+1,LAMBDA+M+JX+1,M-JX);
				AddTimes(X_side,X+JX,LAMBDA+3*M,M-JX);
				AddTimes(X_side,X+JX+1,LAMBDA+3*M,M-JX-1);
				AddTimes(X_side+1,X+JX,LAMBDA+3*M,M-JX);
				Norm(X_side,2.0,M);
				AddTimes(X_side+2*JX,X,LAMBDA+2*JX,M-2*JX);
				AddTimes(X_side,X+2*JX,LAMBDA+4*M,M-2*JX);
				AddTimes(X_side+2*JX,X+1,LAMBDA+2*JX,M-2*JX);
				AddTimes(X_side+2*JX+1,X,LAMBDA+2*JX+1,M-2*JX-1);
				AddTimes(X_side+1,X+2*JX,LAMBDA+4*M+1,M-2*JX);
				AddTimes(X_side,X+2*JX+1,LAMBDA+4*M,M-2*JX-1);
				AddTimes(X_side+JX+2,X,LAMBDA+M+JX+2,M-JX-2);
				AddTimes(X_side+JX,X+2,LAMBDA+M+JX,M-JX);
				AddTimes(X_side,X+JX+2,LAMBDA+3*M,M-JX-2);
				AddTimes(X_side+2,X+JX,LAMBDA+3*M+2,M-JX);
				AddTimes(X_side+2,X,LAMBDA+2*M+2,M-2);
				AddTimes(X_side,X+2,LAMBDA+2*M,M-2);
				Norm(X_side,2.0,M);
				AddTimes(X_side+2*JX+2,X,LAMBDA+2*JX+2,M-2*JX-2);
				AddTimes(X_side,X+2*JX+2,LAMBDA+4*M,M-2*JX-2);
				AddTimes(X_side+2*JX,X+2,LAMBDA+2*JX,M-2*JX);
				AddTimes(X_side+2,X+2*JX,LAMBDA+4*M+2,M-2*JX);
				Norm(X_side,1.0/16.0,M);
			}
		}
	}
}

void LGrad2::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in LGrad2 " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);

	Zero(gs,M); set_bounds(gs_1);
	if (fjc==1) {
		if (geometry=="planar") {
			if (lattice_type=="simple_cubic") { //9 point stencil
				YplusisCtimesX(gs,gs_1,    16.0/36.0,M);
				YplusisCtimesX(gs+1,gs_1,   4.0/36.0,M-1);
				YplusisCtimesX(gs,gs_1+1,   4.0/36.0,M-1);
				YplusisCtimesX(gs+JX,gs_1,  4.0/36.0,M-JX);
				YplusisCtimesX(gs,gs_1+JX,  4.0/36.0,M-JX);
				YplusisCtimesX(gs+JX+1,gs_1,1.0/36.0,M-JX-1);
				YplusisCtimesX(gs+JX,gs_1+1,1.0/36.0,M-JX);
				YplusisCtimesX(gs+1,gs_1+JX,1.0/36.0,M-JX);
				YplusisCtimesX(gs,gs_1+JX+1,1.0/36.0,M-JX-1);
				Times(gs,gs,G1,M);
						 //6point stencil ; classical!
						//Add(gs+JX,gs_1,M-JX);
						//Add(gs,gs_1+JX,M-JX);
						//Add(gs+1,gs_1,M-1);
						//Add(gs,gs_1+1,M-1);
						//Norm(gs,1.0/2.0,M);
						//Add(gs,gs_1,M);
						//Norm(gs,1.0/3.0,M);
						//Times(gs,gs,G1,M);
			} else { //hexagonal //9 point stencil
				YplusisCtimesX(gs,gs_1,    12.0/48.0,M);
				YplusisCtimesX(gs+1,gs_1,   6.0/48.0,M-1);
				YplusisCtimesX(gs,gs_1+1,   6.0/48.0,M-1);
				YplusisCtimesX(gs+JX,gs_1,  6.0/48.0,M-JX);
				YplusisCtimesX(gs,gs_1+JX,  6.0/48.0,M-JX);
				YplusisCtimesX(gs+JX+1,gs_1,3.0/48.0,M-JX-1);
				YplusisCtimesX(gs+JX,gs_1+1,3.0/48.0,M-JX);
				YplusisCtimesX(gs+1,gs_1+JX,3.0/48.0,M-JX);
				YplusisCtimesX(gs,gs_1+JX+1,3.0/48.0,M-JX-1);
				Times(gs,gs,G1,M);
			}
		} else { //geometry==cylindrical. use \lambda's.
			if (lattice_type=="simple cubic") {
				YplusisCtimesX(gs,gs_1,4.0/6.0,M);
				AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
				AddTimes(gs,gs_1+JX,lambda1,M-JX);
				YplusisCtimesX(gs+1,gs_1,1.0/6.0,M-1);
				YplusisCtimesX(gs,gs_1+1,1.0/6.0,M-1);
				Norm(gs,4.0,M);
				AddTimes(gs+JX+1,gs_1,lambda_1+JX+1,M-JX-1);
				AddTimes(gs+JX,gs_1+1,lambda_1+JX,M-JX);
				AddTimes(gs+1,gs_1+JX,lambda1+1,M-JX);
				AddTimes(gs,gs_1+JX+1,lambda1,M-JX-1);
				Norm(gs,1.0/6.0,M);
				Times(gs,gs,G1,M);
					} else {
				YplusisCtimesX(gs,gs_1,2.0/4.0,M);
				AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
				AddTimes(gs,gs_1+JX,lambda1,M-JX);
				YplusisCtimesX(gs+1,gs_1,1.0/4.0,M-1);
				YplusisCtimesX(gs,gs_1+1,1.0/4.0,M-1);
				Norm(gs,2.0,M);
				AddTimes(gs+JX+1,gs_1,lambda_1+JX+1,M-JX-1);
				AddTimes(gs+JX,gs_1+1,lambda_1+JX,M-JX);
				AddTimes(gs+1,gs_1+JX,lambda1+1,M-JX);
				AddTimes(gs,gs_1+JX+1,lambda1,M-JX-1);
				Norm(gs,3.0/12.0,M);
				Times(gs,gs,G1,M);
			}
		}
	}
	if (fjc==2) { //25 point stencil only fjc==2 implemented....
		if (geometry=="planar") { //lattice_type = hexagonal
			Add(gs,gs_1,M);
			Add(gs+JX,gs_1,M-JX);
			Add(gs,gs_1+JX,M-JX);
			Add(gs+1,gs_1,M-1);
			Add(gs,gs_1+1,M-1);
			Add(gs+JX+1,gs_1,M-JX-1);
			Add(gs,gs_1+JX+1,M-JX-1);
			Add(gs+1,gs_1+JX,M-JX);
			Add(gs+JX,gs_1+1,M-JX);
			Norm(gs,2.0,M);
			Add(gs+2*JX,gs_1,M-2*JX);
			Add(gs,gs_1+2*JX,M-2*JX);
			Add(gs+2*JX,gs_1+1,M-2*JX);
			Add(gs+2*JX+1,gs_1,M-2*JX-1);
			Add(gs+1,gs_1+2*JX,M-2*JX);
			Add(gs,gs_1+2*JX+1,M-2*JX-1);
			Add(gs+JX+2,gs_1,M-JX-2);
			Add(gs+JX,gs_1+2,M-JX);
			Add(gs,gs_1+JX+2,M-JX-2);
			Add(gs+2,gs_1+JX,M-JX);
			Add(gs+2,gs_1,M-2);
			Add(gs,gs_1+2,M-2);
			Norm(gs,2.0,M);
			Add(gs+2*JX+2,gs_1,M-2*JX-2);
			Add(gs,gs_1+2*JX+2,M-2*JX-2);
			Add(gs+2*JX,gs_1+2,M-2*JX);
			Add(gs+2,gs_1+2*JX,M-2*JX);
			Norm(gs,1.0/64.0,M);
			Times(gs,gs,G1,M);
		} else { //lattice_type hexagonal.; cylindrical geometry.
			Add(gs,gs_1,M);
			Add(gs+1,gs_1,M-1);
			Add(gs,gs_1+1,M-1);
			Norm(gs,1.0/4.0,M);
			AddTimes(gs+JX,gs_1,LAMBDA+M+JX,M-JX);
			AddTimes(gs+JX+1,gs_1,LAMBDA+M+JX+1,M-JX-1);
			AddTimes(gs+JX,gs_1+1,LAMBDA+M+JX,M-JX);
			AddTimes(gs,gs_1+JX,LAMBDA+3*M,M-JX);
			AddTimes(gs,gs_1+JX+1,LAMBDA+3*M,M-JX-1);
			AddTimes(gs+1,gs_1+JX,LAMBDA+3*M+1,M-JX);
			Norm(gs,2.0,M);
			AddTimes(gs+2*JX,gs_1,LAMBDA+2*JX,M-2*JX);
			AddTimes(gs,gs_1+2*JX,LAMBDA+4*M,M-2*JX);
			AddTimes(gs+2*JX,gs_1+1,LAMBDA+2*JX,M-2*JX);
			AddTimes(gs+2*JX+1,gs_1,LAMBDA+2*JX+1,M-2*JX-1);
			AddTimes(gs+1,gs_1+2*JX,LAMBDA+4*M+1,M-2*JX);
			AddTimes(gs,gs_1+2*JX+1,LAMBDA+4*M,M-2*JX-1);
			AddTimes(gs+JX+2,gs_1,LAMBDA+1*M+JX+2,M-JX-2);
			AddTimes(gs+JX,gs_1+2,LAMBDA+1*M+JX,M-JX);
			AddTimes(gs,gs_1+JX+2,LAMBDA+3*M,M-JX-2);
			AddTimes(gs+2,gs_1+JX,LAMBDA+3*M+2,M-JX);
			AddTimes(gs+2,gs_1,LAMBDA+2*M+2,M-2);
			AddTimes(gs,gs_1+2,LAMBDA+2*M,M-2);
			Norm(gs,2.0,M);
			AddTimes(gs+2*JX+2,gs_1,LAMBDA+2*JX+2,M-2*JX-2);
			AddTimes(gs,gs_1+2*JX+2,LAMBDA+4*M,M-2*JX-2);
			AddTimes(gs+2*JX,gs_1+2,LAMBDA+2*JX,M-2*JX);
			AddTimes(gs+2,gs_1+2*JX,LAMBDA+4*M+2,M-2*JX);
			Norm(gs,1.0/16.0,M);
			Times(gs,gs,G1,M);
		}
	}
}


bool LGrad2::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRange in LGrad2 " << endl;
	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	string diggit;
	bool recognize_keyword;
	//int a; if (range_type=="frozen_range") a=1; else a=0;
	int a=0;
	In[0]->split(range,';',set);
	if (set.size()==2) {
		coor.clear();
		block=true; In[0]->split(set[0],',',coor);
		if (coor.size()!=2) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of two: 'x,y'" << endl; success=false;}
		else {
			diggit=coor[0].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[0]=In[0]->Get_int(coor[0],0); else {
				recognize_keyword=false;
				if (coor[0]=="var_pos") {recognize_keyword=true; r[0]=var_pos; }
				if (coor[0]=="firstlayer") {recognize_keyword=true; r[0] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
				if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[0] < 1-a || r[0] > MX+a) {cout << "In mon " + seg_name + ", for 'pos 1', the x-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MX+a << endl; success =false;}
			diggit=coor[1].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[1]=In[0]->Get_int(coor[1],0); else {
				recognize_keyword=false;
				if (coor[1]=="var_pos") {recognize_keyword=true; r[1]=var_pos;  }
				if (coor[1]=="firstlayer") {recognize_keyword=true; r[1] = 1;}
						//if (coor[1]=="lowerbound") {recognize_keyword=true; r[1] = 0;}
						//if (coor[1]=="upperbound") {recognize_keyword=true; r[1] = MY+1;}
				if (coor[1]=="lastlayer")  {recognize_keyword=true; r[1] = MY;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[1] < 1-a || r[1] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .."<< MY+a << endl; success =false;}
		}
		coor.clear(); In[0]->split(set[1],',',coor);

		if (coor.size()!=2) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of two: 'x,y'" << endl; success=false;}
		else {
			diggit=coor[0].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[3]=In[0]->Get_int(coor[0],0); else {
				recognize_keyword=false;
				if (coor[0]=="var_pos") {recognize_keyword=true; r[3]=var_pos; }
				if (coor[0]=="firstlayer") {recognize_keyword=true; r[3] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
				if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[3] <1-a || r[3] > MX+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the x-coordinate in '" + range_type + "' is out of bounds: " << 1-a << " .. " << MX+a<< endl; success=false;}
			diggit=coor[1].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[4]=In[0]->Get_int(coor[1],0); else {
				recognize_keyword=false;
				if (coor[1]=="var_pos") {recognize_keyword=true; r[4]=var_pos;}
				if (coor[1]=="firstlayer") {recognize_keyword=true; r[4] = 1;}
						//if (coor[1]=="lowerbound") {recognize_keyword=true; r[4] = 0;}
						//if (coor[1]=="upperbound") {recognize_keyword=true; r[4] = MY+1;}
				if (coor[1]=="lastlayer")  {recognize_keyword=true; r[4] = MY;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[4] < 1-a || r[4] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the y-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .."<< MY+a << endl; success =false;}
			if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
			if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
		}
	} else {
		string s;
		In[0]->split(set[0],')',coor);
		s=coor[0].substr(0,1);
		if (s!="(") { //now expect one of the keywords
			block=true;
			cout << "In mon " + seg_name + " and  range_type " + range_type + ", the info was not recognised because when 'gradients>1' the lonely keywords 'firstlayer' 'lastlayers' do not work." << endl;
			success=false;
		} else {
			int px{0},py{0};
			string s;

			if (coor.size()==0)
				block=false;
			else {
				for (size_t i = 0 ; i < coor.size() ; ++i) {
					s=coor[i].substr(1,coor[i].size()-1);
					In[0]->split(s,',',xyz);
					if (xyz.size()!=2) {
						cout << "In mon " + seg_name+ " pinned_range  the expected 'pair of coordinate' -with brackets- structure '(x,y)' was not found. " << endl;  success = false;
					} else {
						px=In[0]->Get_int(xyz[0],0);
						if (px < 1 || px > MX) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}
						py=In[0]->Get_int(xyz[1],0);
						if (py < 1 || py > MY) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 0.." << MY+1 << endl; success =false;}
					}
					H_p[i]=px*JX+fjc-1+py;
				}
			}

		}
	}
	return success;
}

bool LGrad2::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRangeFile in LGrad2 " << endl;
	bool success=true;
	string content;
	vector<string> lines;
	vector<string> sub;
	vector<string> xyz;
	string Infilename=In[0]->name;
	In[0]->split(Infilename,'.',sub);

	int length;
	int length_xyz;
	int px,py,p_i,x,y;
	int i=0;
	if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {
		success=false;
		return success;
	}

	In[0]->split(content,'#',lines);
	length = lines.size();


	if (length == MX*MY) { //expect to read 'mask file';
		i=0;
		if (n_pos==0) {
			while (i<length){
				if (In[0]->Get_int(lines[i],0)==1) n_pos++;
				i++;
			};
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
		} else {
			i=0; p_i=0;
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++)  {
				if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=P(x,y); p_i++;}
				i++;
			}
		}
	} else { //expect to read x,y
		px=0,py=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In[0]->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=2) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'pair of coordinates' 'x,y' was not found. " << endl;  success = false;
				} else {
					px=In[0]->Get_int(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
					py=In[0]->Get_int(xyz[1],0);
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

bool LGrad2::FillMask(int* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
	bool success=true;
	bool readfile=false;
	int length=0;
	int length_px = px.size();

	vector<string> lines;
	int p;
	if (px.size()==0) {
		readfile=true;
		string content;
		success=In[0]->ReadFile(filename,content);
		if (success) {
			In[0]->split(content,'#',lines);
			length = lines.size();
		}
	}
	if (readfile) {
		if (MX*MY!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x,y-directions" << endl;
		} else {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++) Mask[x*JX + y]=In[0]->Get_int(lines[x*JX+y],-1);
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

bool LGrad2::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for LGrad2 " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		for (int x=r[0]; x<r[3]+1; x++)
		for (int y=r[1]; y<r[4]+1; y++)
			H_MASK[P(x,y)]=1;
	} else {
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad2::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	if (geometry !="planar") Dot(result,phi,L,M); 
	else {if (fjc==1) Sum(result,phi,M); else  Dot(result,phi,L,M);}
	return result;
}

void LGrad2::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_bounds(psi);
	Zero(EE,M);
	Real Exmin,Explus,Eymin,Eyplus;
	int x,y,z;
	int r;

	if (geometry=="planar") {
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
				EE[x*MX+y]=pf*(Exmin+Explus+Eymin+Eyplus);
			}
		}
	} else {
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
}


void LGrad2::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y,i;




	Real r;
	Real epsXplus, epsXmin, epsYplus,epsYmin;
	set_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	if (geometry=="planar") {
		C=C*2.0/fjc/fjc;
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				i=x*JX+y;
				epsXmin=eps[i]+eps[i-JX];
				epsXplus=eps[i]+eps[i+JX];
				epsYmin=eps[i]+eps[i-1];
				epsYplus=eps[i]+eps[i+1];
				X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/
				(epsXmin+epsXplus+epsYmin+epsYplus);
			}
		}
	} else {
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
				X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/
				(epsXmin+epsXplus+epsYmin+epsYplus);
			}
		}

	}
	//Cp(psi,X,M);
	YisAminB(g,g,X,M);
   } else { //fixedPsi0 is true
	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++){
			if (Mask[x*JX+y] == 0)
			X[x*JX+y]=0.25*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
			        +psi[x*JX+y-1]  +psi[x*JX+y+1])
				 +0.5*q[x*JX+y]*C/eps[x*JX+y];
		}
	}

	if (geometry=="cylindrical") { //radial geometry in x no radial geometry in y
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++){
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]+=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc;
			}
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


void LGrad2::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
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

	if (geometry=="cylindrical") { //radial geometry in x no radial geometry in y
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++){
				if (Mask[x*JX+y] == 1)
				q[x*JX+y]-=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc*eps[x]/C;
			}
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
if (debug) cout <<"remove_bounds in LGrad2 " << endl;
	int x,y;
	int k;
	if (fjc==1) {
		for (x=0; x<MX+2; x++) {
			X[P(x,0)] = 0;
			X[P(x,MY+1)]=0;
		}
		for (y=0; y<MY+2; y++) {
			X[P(0,y)] = 0;
			X[P(MX+1,y)]=0;
		}
	} else {
		for (x=1-fjc; x<MX+fjc+1; x++) {
			for (k=0; k<fjc; k++) X[P(x,-k)] =0;
			for (k=0; k<fjc; k++) X[P(x,MY+1+k)]=0;
		}
		for (y=1-fjc; y<MY+fjc+1; y++) {
			for (k=0; k<fjc; k++) X[P(-k,y)] = 0;
			for (k=0; k<fjc; k++) X[P(MX+1+k,y)]=0;
		}
	}
}

 
void LGrad2::set_bounds(Real* X){
if (debug) cout <<"set_bounds in LGrad2 " << endl;
	int x,y;
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
		for (x=0; x<1; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
		for (x=MX+1; x<MX+2; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
	} else {
		for (x=1; x<MX+1; x++) {
			for (k=0; k<fjc; k++) {
				X[P(x,-k)] = X[P(x,1+k)];
			}
			for (k=0; k<fjc; k++) {
				X[P(x,MY+1+k)]=X[P(x,MY-k)];
			}
		}
		for (y=1; y<MY+1; y++) {
			for (k=0; k<fjc; k++) {
				X[P(-k,y)] = X[P(1+k,y)];
			}
			for (k=0; k<fjc; k++) {
				X[P(MX+1+k,y)]=X[P(MX-k,y)];
			}
		}
		//corners
		for (x=1-fjc; x<1; x++) {
			for (k=0; k<fjc; k++) {
				X[P(x,-k)] = X[P(x,1+k)];
			}
			for (k=0; k<fjc; k++) {
				X[P(x,MY+1+k)]=X[P(x,MY-k)];
			}
		}
		for (x=MX+1; x<MX+1+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[P(x,-k)] = X[P(x,1+k)];
			}
			for (k=0; k<fjc; k++) {
				X[P(x,MY+1+k)]=X[P(x,MY-k)];
			}
		}
	}
}

void LGrad2::remove_bounds(int *X){
if (debug) cout <<"remove_bounds in LGrad2 " << endl;
	int x,y;
	int k;
	if (fjc==1) {
		for (x=0; x<MX+2; x++) {
			X[P(x,0)] = 0;
			X[P(x,MY+1)]=0;
		}
		for (y=0; y<MY+2; y++) {
			X[P(0,y)] = 0;
			X[P(MX+1,y)]=0;
		}
	} else {
		for (x=1-fjc; x<MX+fjc+1; x++) {
			for (k=0; k<fjc; k++) X[P(x,-k)] =0;
			for (k=0; k<fjc; k++) X[P(x,MY+1+k)]=0;
		}
		for (y=1-fjc; y<MY+fjc+1; y++) {
			for (k=0; k<fjc; k++) X[P(-k,y)] = 0;
			for (k=0; k<fjc; k++) X[P(MX+1+k,y)]=0;
		}
	}
}

 
void LGrad2::set_bounds(int* X){
if (debug) cout <<"set_bounds in LGrad2 " << endl;
	int x,y;
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
		for (x=0; x<1; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
		for (x=MX+1; x<MX+2; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
	} else {
		for (x=1; x<MX+1; x++) {
			for (k=0; k<fjc; k++) {
				X[P(x,-k)] = X[P(x,1+k)];
			}
			for (k=0; k<fjc; k++) {
				X[P(x,MY+1+k)]=X[P(x,MY-k)];
			}
		}
		for (y=1; y<MY+1; y++) {
			for (k=0; k<fjc; k++) {
				X[P(-k,y)] = X[P(1+k,y)];
			}
			for (k=0; k<fjc; k++) {
				X[P(MX+1+k,y)]=X[P(MX-k,y)];
			}
		}
		//corners
		for (x=1-fjc; x<1; x++) {
			for (k=0; k<fjc; k++) {
				X[P(x,-k)] = X[P(x,1+k)];
			}
			for (k=0; k<fjc; k++) {
				X[P(x,MY+1+k)]=X[P(x,MY-k)];
			}
		}
		for (x=MX+1; x<MX+1+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[P(x,-k)] = X[P(x,1+k)];
			}
			for (k=0; k<fjc; k++) {
				X[P(x,MY+1+k)]=X[P(x,MY-k)];
			}
		}
	}
}