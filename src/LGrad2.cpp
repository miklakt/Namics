#include <iostream>
#include <string>
#include "LGrad2.h"
#include "tools.h"

LGrad2::LGrad2(const Input& In_,const std::string& name_): Lattice(In_,name_) {}

LGrad2::~LGrad2() {
NAMICS_DBG("LGrad2 destructor " << std::endl);}

bool LGrad2::CheckLatticeInput(const ParameterStore& parameters) {
	bool success = RejectScalarBoundsInMultiD(parameters);
	success = RejectZBoundsIn2D(parameters) && success;
	success = ReadScaledDimension(parameters, "n_layers_x", MX, 0, "In 'lat' the parameter 'n_layers_x' is required. Problem terminated", "n_layers_x out of bounds, currently: 0.. 1e6; Problem terminated") && success;
	success = ReadScaledDimension(parameters, "n_layers_y", MY, 0, "In 'lat' the parameter 'n_layers_y' is required. Problem terminated", "n_layers_y out of bounds, currently: 0.. 1e6; Problem terminated") && success;

	geometry = parameters.value("geometry", std::string{"planar"});
	success = AssignChoice(geometry, geometry, {"cylindrical", "flat", "planar"}, "In lattice input for 'geometry' not recognized.") && success;
	if (geometry == "flat") geometry = "planar";
	ReadOffsetFirstLayer(parameters);
	success = ReadBoundaryCondition(parameters, "lowerbound_x", 0, {"mirror", "surface"}, "for 'lowerbound_x' boundary condition not recognized.  ") && success;
	success = ReadBoundaryCondition(parameters, "upperbound_x", 3, {"mirror", "surface"}, "for 'upperbound_x' boundary condition not recognized. ") && success;
	success = ReadBoundaryCondition(parameters, "lowerbound_y", 1, {"mirror", "surface"}, "for 'lowerbound_y' boundary condition not recognized. ") && success;
	success = ReadBoundaryCondition(parameters, "upperbound_y", 4, {"mirror", "surface"}, "for 'upperbound_Y' boundary condition not recognized. ") && success;
	return success;
}


void LGrad2:: ComputeLambdas() {
	Real r, VL, LS;
	Real rlow, rhigh;


	//	}
		if (fjc==1) {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++) {
				r=offset_first_layer + 1.0*x;
				lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]*lambda;
				lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]*lambda;
				lambda0[P(x,y)]=1.0-2.0*lambda;
				if (fcc_sites) {
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
}

bool LGrad2::PutM() {
NAMICS_DBG("PutM in LGrad2 " << std::endl);	bool success=true;

	if (geometry=="cylindrical")
		volume = MY*PIE*(std::pow(MX+offset_first_layer,2)-std::pow(offset_first_layer,2));
	else volume = MX*MY;
	JX=MY+2*fjc; JY=1; JZ=0; M=(MX+2*fjc)*(MY+2*fjc);

	Accesible_volume=volume;
	return success;
}

Real LGrad2:: Moment(Real* X,Real Xb, int n) {
NAMICS_DBG("Moment in LGrad2 " << std::endl);	Real Result=0;
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
			if (zz>0) Result+= std::pow(zz,n)*Nz;
		}
	} else {
	}
	return Result/fjc;
}

Real LGrad2::WeightedSum(Real* X){
NAMICS_DBG("weighted sum in LGrad2 " << std::endl);	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		sum = std::accumulate(X, X + M, Real{0}) / (fjc * fjc);
	} else	{
		sum = std::inner_product(X, X + M, L.begin(), Real{0});
	}
	return sum;
}

void LGrad2::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG(" Side in LGrad2 " << std::endl);	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	Real* fcc_lambda_1 = this->fcc_lambda_1.data();
	Real* fcc_lambda1 = this->fcc_lambda1.data();
	Real* lambda_1 = this->lambda_1.data();
	Real* lambda1 = this->lambda1.data();
	Real* LAMBDA = this->LAMBDA.data();
	std::fill_n(X_side, M, 0);//set_bounds(X);

	if (fcc_sites) {
		Real C1 = 1.0 / 3.0;
		add_shifted(X_side, X, M, C1);
		add_weighted(X_side + JX, X, fcc_lambda_1 + JX, M - JX);
		add_weighted(X_side, X + JX, fcc_lambda1, M - JX);
		add_shifted(X_side + 1, X, M - 1, C1);
		add_shifted(X_side, X + 1, M - 1, C1);
		add_weighted(X_side + JX + 1, X, fcc_lambda_1 + JX + 1, M - JX - 1);
		add_weighted(X_side + JX, X + 1, fcc_lambda_1 + JX, M - JX - 1);
		add_weighted(X_side + 1, X + JX, fcc_lambda1 + 1, M - JX - 1);
		add_weighted(X_side, X + JX + 1, fcc_lambda1, M - JX - 1);
		scale_span(X_side, M, C1);

	} else {
		if (fjc==1) {
			if (lattice_type ==simple_cubic) {
				Real C1 = 4.0 / 6.0;
				Real C2 = 1.0 / 6.0;
				Real C3 = 4.0;
				add_shifted(X_side, X, M, C1);
				add_weighted(X_side + JX, X, lambda_1 + JX, M - JX);
				add_weighted(X_side, X + JX, lambda1, M - JX);
				add_shifted(X_side + 1, X, M - 1, C2);
				add_shifted(X_side, X + 1, M - 1, C2);
				scale_span(X_side, M, C3);
				add_weighted(X_side + JX + 1, X, lambda_1 + JX + 1, M - JX - 1);
				add_weighted(X_side + JX, X + 1, lambda_1 + JX, M - JX - 1);
				add_weighted(X_side + 1, X + JX, lambda1 + 1, M - JX - 1);
				add_weighted(X_side, X + JX + 1, lambda1, M - JX - 1);
				scale_span(X_side, M, C2);
			} else {
				Real C1 = 2.0 / 4.0;
				Real C2 = 1.0 / 4.0;
				Real C3 = 2.0;
				add_shifted(X_side, X, M, C1);
				add_weighted(X_side + JX, X, lambda_1 + JX, M - JX);
				add_weighted(X_side, X + JX, lambda1, M - JX);
				add_shifted(X_side + 1, X, M - 1, C2);
				add_shifted(X_side, X + 1, M - 1, C2);
				scale_span(X_side, M, C3);
				add_weighted(X_side + JX + 1, X, lambda_1 + JX + 1, M - JX - 1);
				add_weighted(X_side + JX, X + 1, lambda_1 + JX, M - JX - 1);
				add_weighted(X_side + 1, X + JX, lambda1 + 1, M - JX - 1);
				add_weighted(X_side, X + JX + 1, lambda1, M - JX - 1);
				scale_span(X_side, M, C2);
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
				if (block !=1) scale_span(X_side, M, 2.0); else scale_span(X_side, M, C);
			}
		}
	}
}


void LGrad2::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
NAMICS_DBG(" propagate in LGrad2 " << std::endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	Real* lambda_1 = this->lambda_1.data();
	Real* lambda1 = this->lambda1.data();
	Real* LAMBDA = this->LAMBDA.data();
	std::fill_n(gs, M, 0); set_bounds(gs_1);
	if (fjc==1) {
		if (lattice_type==simple_cubic) {
			Real C1=4.0/6.0;
			Real C2=1.0/6.0;
			Real C3=4.0;
			add_shifted(gs, gs_1, M, C1);
			add_weighted(gs + JX, gs_1, lambda_1 + JX, M - JX);
			add_weighted(gs, gs_1 + JX, lambda1, M - JX);
			add_shifted(gs + 1, gs_1, M - 1, C2);
			add_shifted(gs, gs_1 + 1, M - 1, C2);
			scale_span(gs, M, C3);
			add_weighted(gs + JX + 1, gs_1, lambda_1 + JX + 1, M - JX - 1);
			add_weighted(gs + JX, gs_1 + 1, lambda_1 + JX, M - JX - 1);
			add_weighted(gs + 1, gs_1 + JX, lambda1 + 1, M - JX - 1);
			add_weighted(gs, gs_1 + JX + 1, lambda1, M - JX - 1);
			scale_span(gs, M, C2);
			std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });
		} else { //9 point stencil; hexagonal
			Real C1=0.5;
			Real C2=0.25;
			Real C3=2.0;
			add_shifted(gs, gs_1, M, C1);
			add_weighted(gs + JX, gs_1, lambda_1 + JX, M - JX);
			add_weighted(gs, gs_1 + JX, lambda1, M - JX);
			add_shifted(gs + 1, gs_1, M - 1, C2);
			add_shifted(gs, gs_1 + 1, M - 1, C2);
			scale_span(gs, M, C3);
			add_weighted(gs + JX + 1, gs_1, lambda_1 + JX + 1, M - JX - 1);
			add_weighted(gs + JX, gs_1 + 1, lambda_1 + JX, M - JX - 1);
			add_weighted(gs + 1, gs_1 + JX, lambda1 + 1, M - JX - 1);
			add_weighted(gs, gs_1 + JX + 1, lambda1, M - JX - 1);
			scale_span(gs, M, C2);
			std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });
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
		std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });


	}
}


void LGrad2::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	if (geometry == "planar") {
		set_M_bounds(psi);
		std::fill_n(EE, M, 0);
		Real Exmin,Explus,Eymin,Eyplus;
		int x,y,z;

		pf = pf / 2.0;
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
		return;
	}
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
	if (geometry == "planar") {
		if (!fixedPsi0) {
			C = C * 2.0 / fjc / fjc;
			for (x=fjc; x<MX+fjc; x++) {
				for (y=fjc; y<MY+fjc; y++) {
					i=x*JX+y;
					Real epsXmin=eps[i]+eps[i-JX];
					Real epsXplus=eps[i]+eps[i+JX];
					Real epsYmin=eps[i]+eps[i-1];
					Real epsYplus=eps[i]+eps[i+1];
					if (x==fjc) a=psi[i-JX]; else a=X[i-JX];
					if (y==fjc) b=psi[i-1]; else b=X[i-1];
					X[i]= (C*q[i]+epsXmin*a+epsXplus*psi[i+JX]+epsYmin*b+epsYplus*psi[i+1])/
					(epsXmin+epsXplus+epsYmin+epsYplus);
				}
			}
			std::transform(g, g + M, X.data(), g, [](auto a, auto b) { return a - b; });
		} else {
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
		return;
	}
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
		std::transform(g, g + M, X.data(), g, [](auto a, auto b) { return a - b; });
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
	if (geometry == "planar") {
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
		return;
	}
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
NAMICS_DBG("remove_bounds in LGrad2 " << std::endl);	int x,y;
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
NAMICS_DBG("set_bounds_x XY in LGrad2 " << std::endl);
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
			std::cout <<"set_bounds_x error" << std::endl;
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
NAMICS_DBG("set_bounds_y XY in LGrad2 " << std::endl);
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
NAMICS_DBG("set_bounds_x X in LGrad2 " << std::endl);	int y;
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
NAMICS_DBG("set_bounds_y X in LGrad2 " << std::endl);	int x;
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
NAMICS_DBG("set_bounds in LGrad2 " << std::endl);	int x,y;
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
				X[(MX+1)*JX+      0] = X[(MX+1)*JX+BY1];
				X[(MX+1)*JX+   MY+1]  =X[(MX+1)*JX+BYM];
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
		}
	}

void LGrad2::set_M_bounds(Real* X){
NAMICS_DBG("set_bounds in LGrad2 " << std::endl);	int x,y;
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
NAMICS_DBG("remove_bounds in LGrad2 " << std::endl);	int x,y;
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
NAMICS_DBG("set_bounds in LGrad2 " << std::endl);	int x,y;
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
				X[(MX+1)*JX+0] = X[(MX+1)*JX+BY1];
				X[(MX+1)*JX+MY+1]=X[(MX+1)*JX+BYM];
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

Real LGrad2::ComputeGN(Real* G, int M){
	return WeightedSum(G);
}

void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb){
	for (int __i = 0; __i < M; ++__i) (phi)[__i] += (Gf)[__i] * (Gb)[__i];
}

void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb,Real degeneracy){
	for (int __i = 0; __i < M; ++__i) (phi)[__i] += degeneracy * (Gf)[__i] * (Gb)[__i];
}

void LGrad2::Initiate(Real* G,Real* Gz){
	std::copy_n(Gz, M, G);
}

void LGrad2::Terminate(Real* Gz,Real* G){
	std::copy_n(G, M, Gz);
}

bool LGrad2:: PutMask(Real* MASK,std::vector<int>px,std::vector<int>py,std::vector<int>pz,int R){
	(void)pz;
NAMICS_DBG("PutMask in LGrad2 " << std::endl);	//R*=fjc; //is already done in segment
	if (geometry == "planar") {
		(void)R;
		(void)py;
		(void)px;
		(void)MASK;
		bool success=false;
		std::cout <<"PutMask does not make sense in planar 2 gradient system " << std::endl;
		return success;
	}
	bool success=true;
	int length =px.size();
	int X,Y;
	int dx,dy;
	Real teller,noemer;
	if (length > 1) {
		std::cout <<"In two gradient system, we can have just one particle: we found " <<length <<"particles. " << std::endl;
		return false;
	}
	for (int i =0; i<length; i++) {
		int xx,yy;
		xx=px[i]; yy=py[i];
		if (xx !=0) {
			std::cout <<"In two gradients system, we expect the particle at the central axis" << std::endl;
			return false;
		}
		if (R>MX || R>MY) {std::cout <<" particle should be smaller than size of box in X or Y direction" << std::endl; return false;}
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
