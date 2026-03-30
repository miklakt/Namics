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


void LGrad2::LReflect(Real *H, Real *P, Real *Q) {
	Real* l_1 = this->l_1.data();
	Real* l_11 = this->l_11.data();
	for (int __i = 0; __i < (M-JX); ++__i) (H)[__i] = (l_1+JX)[__i] * (P)[__i];
	for (int __i = 0; __i < (M-JX); ++__i) (H)[__i] += (l_11+JX)[__i] * (Q+JX)[__i];
}

void LGrad2::UReflect(Real *H, Real *P, Real *Q) {
	Real* l1 = this->l1.data();
	Real* l11 = this->l11.data();
	for (int __i = 0; __i < (M-JX); ++__i) (H+JX)[__i] = (l1)[__i] * (P+JX)[__i];
	for (int __i = 0; __i < (M-JX); ++__i) (H+JX)[__i] += (l11)[__i] * (Q)[__i];
}


void LGrad2::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	Real* H = this->H.data();
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
			add_terms(gx0 + JX, M - JX, {{gz3, P[1]}, {gz4, P[1]}, {gz5, P[1]}, {gz6, P[1]}, {gz7, P[1]}, {gz8, P[1]}});
			add_terms(gx11, M - JX, {{gz3 + JX, P[1]}, {gz4 + JX, P[1]}, {gz5 + JX, P[1]}, {gz6 + JX, P[1]}, {gz7 + JX, P[1]}, {gz8 + JX, P[1]}});


			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			add_terms(gx3 + JY, M - JY, {{gz0, P[1]}, {gz1, P[1]}, {gz2, P[1]}, {gz3, P[0]}, {gz4, P[0]}, {gz5, P[0]}, {gz9, P[1]}, {gz10, P[1]}, {gz11, P[1]}});
			add_terms(gx8, M - JY, {{gz0 + JY, P[1]}, {gz1 + JY, P[1]}, {gz2 + JY, P[1]}, {gz6 + JY, P[0]}, {gz7 + JY, P[0]}, {gz8 + JY, P[0]}, {gz9 + JY, P[1]}, {gz10 + JY, P[1]}, {gz11 + JY, P[1]}});

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			add_terms(gx5, M, {{gz0, P[1]}, {gz1, P[1]}, {gz2, P[1]}, {gz3, P[0]}, {gz4, P[0]}, {gz5, P[0]}, {gz9, P[1]}, {gz10, P[1]}, {gz11, P[1]}});
			add_terms(gx6, M, {{gz0, P[1]}, {gz1, P[1]}, {gz2, P[1]}, {gz6, P[0]}, {gz7, P[0]}, {gz8, P[0]}, {gz9, P[1]}, {gz10, P[1]}, {gz11, P[1]}});

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			add_terms(gx1 + JX, M - JX, {{gz3, P[1]}, {gz4, P[1]}, {gz5, P[1]}, {gz6, P[1]}, {gz7, P[1]}, {gz8, P[1]}});
			add_terms(gx10, M - JX, {{gz3 + JX, P[1]}, {gz4 + JX, P[1]}, {gz5 + JX, P[1]}, {gz6 + JX, P[1]}, {gz7 + JX, P[1]}, {gz8 + JX, P[1]}});

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz0,gz11,0); set_bounds_y(gz1,gz10,0); set_bounds_y(gz5,gz6,0);

			add_terms(gx2, M - JY, {{gz0 + JY, P[0]}, {gz1 + JY, P[0]}, {gz2 + JY, P[0]}, {gz3 + JY, P[1]}, {gz4 + JY, P[1]}, {gz5 + JY, P[1]}, {gz6 + JY, P[1]}, {gz7 + JY, P[1]}, {gz8 + JY, P[1]}});
			add_terms(gx9 + JY, M - JY, {{gz3, P[1]}, {gz4, P[1]}, {gz5, P[1]}, {gz6, P[1]}, {gz7, P[1]}, {gz8, P[1]}, {gz9, P[0]}, {gz10, P[0]}, {gz11, P[0]}});

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			add_terms(gx4 + JY, M - JY, {{gz0, P[1]}, {gz1, P[1]}, {gz2, P[1]}, {gz3, P[0]}, {gz4, P[0]}, {gz5, P[0]}, {gz9, P[1]}, {gz10, P[1]}, {gz11, P[1]}});
			add_terms(gx7, M - JY, {{gz0 + JY, P[1]}, {gz1 + JY, P[1]}, {gz2 + JY, P[1]}, {gz6 + JY, P[0]}, {gz7 + JY, P[0]}, {gz8 + JY, P[0]}, {gz9 + JY, P[1]}, {gz10 + JY, P[1]}, {gz11 + JY, P[1]}});

			for (int k=0; k<12; k++) std::transform(gs+k*M, gs+(k+1)*M, g, gs+k*M, [](auto a, auto b) { return a * b; });



		} else {//simple _cubic should work
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M,*gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M,*gx4=gs+4*M;
			Real *g=G1;

			std::fill_n(gs, 5*M, 0);
			LReflect(H,gz0,gz4);
			add_terms(gx0 + JX, M - JX, {{H, P[0]}, {gz1, P[1]}, {gz2, 2 * P[1]}, {gz3, P[1]}});
			add_terms(gx1 + JY, M - JY, {{gz0, P[1]}, {gz1, P[0]}, {gz2, 2 * P[1]}, {gz4, P[1]}});
			add_terms(gx2, M, {{gz0, P[1]}, {gz1, P[1]}, {gz2, P[0]}, {gz3, P[1]}, {gz4, P[1]}});
			add_terms(gx3, M - JY, {{gz0 + JY, P[1]}, {gz2 + JY, 2 * P[1]}, {gz3 + JY, P[0]}, {gz4 + JY, P[1]}});
			add_terms(gx4, M - JX, {{gz1 + JX, P[1]}, {gz2 + JX, 2 * P[1]}, {gz3 + JX, P[1]}});
			UReflect(H,gz4,gz0);
			add_from_source(H + JX, M - JX, {{gx4, P[0]}});

			for (int k=0; k<5; k++) std::transform(gs+k*M, gs+(k+1)*M, g, gs+k*M, [](auto a, auto b) { return a * b; });
		}
	} else {
		if (lattice_type==hexagonal) {
			std::cout <<"stencil_full, hexagonal, Markov 2, cyl coordinates not implemented" << std::endl;
		} else {
			std::cout <<"stencil_full, simple_cubic, Markov 2, cyl coordinates not implemented" << std::endl;
		}
	}
}

void LGrad2::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	Real* H = this->H.data();

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
			add_from_source(H, M - JX, {{gx3 + JX, P[1]}, {gx4 + JX, P[1]}, {gx5 + JX, P[1]}, {gx6 + JX, P[1]}, {gx7 + JX, P[1]}, {gx8 + JX, P[1]}, {gx9 + JX, P[0]}, {gx10 + JX, P[0]}, {gx11 + JX, P[0]}});

			UReflect(H,gz0,gz11);
			add_from_source(H + JX, M - JX, {{gx0, P[0]}, {gx1, P[0]}, {gx2, P[0]}, {gx3, P[1]}, {gx4, P[1]}, {gx5, P[1]}, {gx6, P[1]}, {gx7, P[1]}, {gx8, P[1]}});

			remove_bounds(gz0);remove_bounds(gz11);
			set_bounds_y(gz3,gz8,0);

			add_from_source(gz8, M - JY, {{gx0 + JY, P[1]}, {gx1 + JY, P[1]}, {gx2 + JY, P[1]}, {gx6 + JY, P[0]}, {gx7 + JY, P[0]}, {gx8 + JY, P[0]}, {gx9 + JY, P[1]}, {gx10 + JY, P[1]}, {gx11 + JY, P[1]}});
			add_from_source(gz3 + JY, M - JY, {{gx0, P[1]}, {gx1, P[1]}, {gx2, P[1]}, {gx3, P[0]}, {gx4, P[0]}, {gx5, P[0]}, {gx9, P[1]}, {gx10, P[1]}, {gx11, P[1]}});

			remove_bounds(gz3);remove_bounds(gz8);

			add_from_source(gz6, M, {{gx0, P[1]}, {gx1, P[1]}, {gx2, P[1]}, {gx6, P[0]}, {gx7, P[0]}, {gx8, P[0]}, {gx9, P[1]}, {gx10, P[1]}, {gx11, P[1]}});
			add_from_source(gz5, M, {{gx0, P[1]}, {gx1, P[1]}, {gx2, P[1]}, {gx3, P[0]}, {gx4, P[0]}, {gx5, P[0]}, {gx9, P[1]}, {gx10, P[1]}, {gx11, P[1]}});

			remove_bounds(gz5); remove_bounds(gz6);
			set_bounds_x(gz1,gz10,0);

			LReflect(H,gz10,gz1);

			add_from_source(H, M - JX, {{gx3 + JX, P[1]}, {gx4 + JX, P[1]}, {gx5 + JX, P[1]}, {gx6 + JX, P[1]}, {gx7 + JX, P[1]}, {gx8 + JX, P[1]}, {gx9 + JX, P[0]}, {gx10 + JX, P[0]}, {gx11 + JX, P[0]}});

			UReflect(H,gz1,gz10);

			add_from_source(H + JX, M - JX, {{gx0, P[0]}, {gx1, P[0]}, {gx2, P[0]}, {gx3, P[1]}, {gx4, P[1]}, {gx5, P[1]}, {gx6, P[1]}, {gx7, P[1]}, {gx8, P[1]}});

			remove_bounds(gz1);remove_bounds(gz10);
			add_from_source(H + JY, M - JY, {{gx3, P[1]}, {gx4, P[1]}, {gx5, P[1]}, {gx6, P[1]}, {gx7, P[1]}, {gx8, P[1]}, {gx9, P[0]}, {gx10, P[0]}, {gx11, P[0]}});
			add_from_source(H, M - JY, {{gx0 + JY, P[0]}, {gx1 + JY, P[0]}, {gx2 + JY, P[0]}, {gx3 + JY, P[1]}, {gx4 + JY, P[1]}, {gx5 + JY, P[1]}, {gx6 + JY, P[1]}, {gx7 + JY, P[1]}, {gx8 + JY, P[1]}});

			remove_bounds(gz2);remove_bounds(gz9);

			add_from_source(gz7, M - JY, {{gx0 + JY, P[1]}, {gx1 + JY, P[1]}, {gx2 + JY, P[1]}, {gx6 + JY, P[0]}, {gx7 + JY, P[0]}, {gx8 + JY, P[0]}, {gx9 + JY, P[1]}, {gx10 + JY, P[1]}, {gx11 + JY, P[1]}});
			add_from_source(gz4 + JY, M - JY, {{gx0, P[1]}, {gx1, P[1]}, {gx2, P[1]}, {gx3, P[0]}, {gx4, P[0]}, {gx5, P[0]}, {gx9, P[1]}, {gx10, P[1]}, {gx11, P[1]}});

			for (int k=0; k<12; k++) std::transform(gs+k*M, gs+(k+1)*M, g, gs+k*M, [](auto a, auto b) { return a * b; });


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
			add_from_source(H, M - JX, {{gx1 + JX, P[1]}, {gx2 + JX, P[1]}, {gx3 + JX, P[1]}, {gx4 + JX, P[0]}});

			add_from_source(gz3, M - JY, {{gx0 + JY, P[1]}, {gx2 + JY, P[1]}, {gx3 + JY, P[0]}, {gx4 + JY, P[1]}});

			add_from_source(gz2, M, {{gx0, 2 * P[1]}, {gx1, 2 * P[1]}, {gx2, P[0]}, {gx3, 2 * P[1]}, {gx4, 2 * P[1]}});

			add_from_source(gz1 + JY, M - JY, {{gx0, P[1]}, {gx1, P[0]}, {gx2, P[1]}, {gx4, P[1]}});

			UReflect(H,gz0,gz4);
			add_from_source(H + JX, M - JX, {{gx0, P[0]}, {gx1, P[1]}, {gx2, P[1]}, {gx3, P[1]}});

			for (int k=0; k<5; k++) std::transform(gs+k*M, gs+(k+1)*M, g, gs+k*M, [](auto a, auto b) { return a * b; });
		}
	} else {
		if (lattice_type==hexagonal) {
			std::cout <<"stencil_full, cyl coordinates, hexagonal, Markov 2 not implemented" << std::endl;
		} else {
			std::cout <<"stencil_full, cyl coordinates, simple_cubic, Markov 2 not implemented" << std::endl;
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

Real LGrad2::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	if (Markov==2) {
		if (lattice_type == hexagonal && !stencil_full) {
			for (int k=0; k<12; k++) {
				GN += WeightedSum(G+k*M);
			}
			GN /=12.0;
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
NAMICS_DBG("LGrad2:: terminate " << std::endl);	if (Markov==2) {
		std::cout <<"terminate in markov==2 is not tested" << std::endl;
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
		std::cout <<"possible problem in LGrad2::Terminate " << std::endl;
	} else std::copy_n(G, M, Gz);
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
