#include "LGrad1.h"
#include "tools.h"

LGrad1::LGrad1(const Input& In_,const std::string& name_): Lattice(In_,name_) {
NAMICS_DBG("LGrad1 constructor " << std::endl);}

bool LGrad1::CheckLatticeInput(const ParameterStore& parameters) {
	bool success = ReadScaledDimension(parameters, "n_layers", MX, 0, "In 'lat' the parameter 'n_layers' is required. Problem terminated", "n_layers out of bounds, currently: 0..1e6; Problem terminated");
	success = RejectAxisBoundsIn1D(parameters) && success;

	geometry = parameters.value("geometry", std::string{"planar"});
	success = AssignChoice(geometry, geometry, {"spherical", "cylindrical", "flat", "planar"}, "In lattice input for 'geometry' not recognized.") && success;
	if (geometry == "flat") geometry = "planar";
	ReadOffsetFirstLayer(parameters);
	success = ReadBoundaryCondition(parameters, "lowerbound", 0, {"mirror", "surface", "periodic"}, "For 'lowerbound' boundary condition not recognized. ") && success;
	success = ReadBoundaryCondition(parameters, "upperbound", 3, {"mirror", "surface", "periodic"}, "For 'upperbound' boundary condition not recognized.") && success;
	return success;
}

void LGrad1:: ComputeLambdas() {
NAMICS_DBG("LGrad1 computeLambda's " << std::endl);
	Real r, VL, LS;
	Real rlow, rhigh;

	if (geometry == "planar") {
		std::fill_n(L.begin() + 1, MX, Real{1});
		if (fjc == 1) {
			std::fill_n(lambda1.begin() + 1, MX, lambda);
			std::fill_n(lambda_1.begin() + 1, MX, lambda);
			std::fill_n(lambda0.begin() + 1, MX, 1.0 - 2.0 * lambda);
		} else {
			for (int i = 0; i < M; ++i) {
				L[i] = 1.0 / fjc;
				LAMBDA[i] = 1.0 / (2 * (FJC - 1));
				LAMBDA[i + (FJC - 1) * M] = 1.0 / (2 * (FJC - 1));
				LAMBDA[i + (FJC - 1) / 2 * M] = 1.0 / (FJC - 1);
				for (int j = 1; j < FJC / 2; j++) {
					LAMBDA[i + j * M] = 1.0 / (FJC - 1);
					LAMBDA[i + (FJC - j - 1) * M] = 1.0 / (FJC - 1);
				}
			}
		}
		return;
	}

	if (fcc_sites){
		if (geometry=="cylindrical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=PIE*(std::pow(r,2)-std::pow(r-1,2));
				lambda1[i]=2.0*PIE*r/L[i]/3.0;
				lambda_1[i]=2.0*PIE*(r-1)/L[i]/3.0;
				lambda0[i]=1.0/3.0;
			}
		}
		if (geometry=="spherical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=4.0/3.0*PIE*(std::pow(r,3)-std::pow(r-1,3));
				fcc_lambda1[i]=4.0*PIE*std::pow(r,2)/L[i]/3.0;
				fcc_lambda_1[i]=4.0*PIE*std::pow(r-1,2)/L[i]/3.0;
				fcc_lambda0[i]=1.0-fcc_lambda1[i]-fcc_lambda_1[i];
			}
		}
	}

	if (fjc==1) {
		if (geometry=="cylindrical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=PIE*(std::pow(r,2)-std::pow(r-1,2));
				lambda1[i]=2.0*PIE*r/L[i]*lambda;
				lambda_1[i]=2.0*PIE*(r-1)/L[i]*lambda;
				lambda0[i]=1.0-2.0*lambda;
			}
		}
		if (geometry=="spherical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=4.0/3.0*PIE*(std::pow(r,3)-std::pow(r-1,3));
				lambda1[i]=4.0*PIE*std::pow(r,2)/L[i]*lambda;
				lambda_1[i]=4.0*PIE*std::pow(r-1,2)/L[i]*lambda;
				lambda0[i]=1.0-lambda1[i]-lambda_1[i];

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
	}
}

void LGrad1::PutM() {
NAMICS_DBG("PutM in LGrad1 " << std::endl);
	JX=1; JY=0; JZ=0; M=MX+2*fjc;
	if (geometry=="planar") {volume = MX/fjc; }
	if (geometry=="spherical") {volume = 4.0/3.0*PIE*(std::pow(MX+offset_first_layer,3)-std::pow(offset_first_layer,3))/fjc/fjc/fjc;}
	if (geometry=="cylindrical") {volume = PIE*(std::pow(MX+offset_first_layer,2)-std::pow(offset_first_layer,2))/fjc/fjc;}

	Accesible_volume=volume;
}

Real LGrad1:: Moment(Real* X,Real Xb, int n) {
NAMICS_DBG("Moment in LGrad1 " << std::endl);	Real Result=0;
	Real cor;
	remove_bounds(X);
	for (int i = fjc; i<M; i++) {
		cor = (i-fjc+0.5)/fjc;
		Result += std::pow(cor,n)*(X[i]-Xb)*L[i];
	}
	return Result/fjc;
}

Real LGrad1::WeightedSum(Real* X){
NAMICS_DBG("weighted sum in LGrad1 " << std::endl);	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		sum = std::accumulate(X, X + M, Real{0}) / fjc;
	} else {
		sum = std::inner_product(X, X + M, L.begin(), Real{0});
	}
	return sum;
}

void LGrad1::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG(" Side in LGrad1 " << std::endl);
	Real* fcc_lambda0 = this->fcc_lambda0.data();
	Real* fcc_lambda_1 = this->fcc_lambda_1.data();
	Real* fcc_lambda1 = this->fcc_lambda1.data();
	Real* lambda0 = this->lambda0.data();
	Real* lambda_1 = this->lambda_1.data();
	Real* lambda1 = this->lambda1.data();
	Real* LAMBDA = this->LAMBDA.data();
	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0);//set_bounds(X);
	int j, kk;

	if (fcc_sites) {
		add_weighted(X_side, X, fcc_lambda0, M);
		add_weighted(X_side + 1, X, fcc_lambda_1 + 1, M - 1);
		add_weighted(X_side, X + 1, fcc_lambda1, M - 1);

	} else {
		if (fjc==1) {
			add_weighted(X_side, X, lambda0, M);
			add_weighted(X_side + 1, X, lambda_1 + 1, M - 1);
			add_weighted(X_side, X + 1, lambda1, M - 1);
		} else {

			for (j = 0; j < FJC/2; j++) {
				kk = (FJC-1)/2-j;
				add_weighted(X_side + kk, X, LAMBDA + j * M + kk, M - kk);
				add_weighted(X_side, X + kk, LAMBDA + (FJC - j - 1) * M, M - kk);
			}
			add_weighted(X_side, X, LAMBDA + (FJC - 1) / 2 * M, M);

		}
	}
}



void LGrad1::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
NAMICS_DBG(" propagate in LGrad1 " << std::endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	Real* lambda0 = this->lambda0.data();
	Real* lambda_1 = this->lambda_1.data();
	Real* lambda1 = this->lambda1.data();
	Real* LAMBDA = this->LAMBDA.data();
	int kk;
	int j;
	std::fill_n(gs, M, 0); set_bounds(gs_1);

	if (fjc==1) {
		add_weighted(gs, gs_1, lambda0, M);
		add_weighted(gs + 1, gs_1, lambda_1 + 1, M - 1);
		add_weighted(gs, gs_1 + 1, lambda1, M - 1);
		std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });

	} else {
		for (j = 0; j < FJC/2; j++) {
			kk = (FJC-1)/2-j;
			add_weighted(gs + kk, gs_1, LAMBDA + j * M + kk, M - kk);
			add_weighted(gs, gs_1 + kk, LAMBDA + (FJC - j - 1) * M, M - kk);
		}
		add_weighted(gs, gs_1, LAMBDA + (FJC - 1) / 2 * M, M);
		std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });
	}
}


void LGrad1::UpdateEE(Real* EE, Real* psi) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	if (geometry == "planar") {
		set_M_bounds(psi);
		std::fill_n(EE, M, 0);
		Real Exmin, Explus;
		pf = pf / 2.0 * fjc * fjc;
		Explus = psi[fjc - 1] - psi[fjc];
		Explus *= Explus;
		for (int x = fjc; x < MX + fjc; x++) {
			Exmin = Explus;
			Explus = psi[x] - psi[x + 1];
			Explus *= Explus;
			EE[x] = pf * (Exmin + Explus);
		}
		return;
	}
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
	Real C =e*e/(eps0*k_BT*bond_length);

	if (geometry == "planar") {
		if (!fixedPsi0) {
			C = C * 2.0 / fjc / fjc;
			epsXplus = eps[fjc - 1] + eps[fjc];
			a = 0; b = psi[fjc - 1]; c = psi[fjc];
			for (x = fjc; x < MX + fjc; x++) {
				epsXmin = epsXplus;
				epsXplus = eps[x] + eps[x + 1];
				if (x == fjc) a = psi[fjc - 1]; else a = X[x - 1];
				X[x] = (epsXmin * a + C * q[x] + epsXplus * psi[x + 1]) / (epsXmin + epsXplus);
			}
			std::transform(g, g + M, X.data(), g, [](auto a, auto b) { return a - b; });
		} else {
			a = 0; b = psi[fjc - 1]; c = psi[fjc];
			for (x = fjc; x < MX + fjc; x++) {
				a = b; b = c; c = psi[x + 1];
				if (Mask[x] == 0) psi[x] = 0.5 * (a + c) + q[x] * C / eps[x];
			}
			if (grad_epsilon) {
				a = 0; b = psi[fjc - 1]; c = psi[fjc]; a_ = 0; b_ = eps[fjc - 1]; c_ = eps[fjc];
				for (x = fjc; x < MX + fjc; x++) {
					a = b; b = c; c = psi[x + 1]; a_ = b_; b_ = c_; c_ = eps[x + 1];
					if (Mask[x] == 0) {
						psi[x] += 0.25 * (c_ - a_) * (c - a) / eps[x] * fjc * fjc;
					}
				}
			}
			for (x = fjc; x < MX + fjc; x++)
			if (Mask[x] == 0) {
				g[x] -= psi[x];
			}
		}
		return;
	}

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
		 }
	}
	std::transform(g, g + M, X.data(), g, [](auto a, auto b) { return a - b; });
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
	if (geometry == "planar") {
		a = 0; b = psi[fjc - 1]; c = psi[fjc];
		for (x = fjc; x < MX + fjc; x++) {
			a = b; b = c; c = psi[x + 1];
			if (Mask[x] == 1) q[x] = -0.5 * (a - 2 * b + c) * fjc * fjc * eps[x] / C;
		}
		if (grad_epsilon) {
			a = 0; b = psi[fjc - 1]; c = psi[fjc]; a_ = 0; b_ = eps[fjc - 1]; c_ = eps[fjc];
			for (x = fjc; x < MX + fjc; x++) {
				a = b; b = c; c = psi[x + 1]; a_ = b_; b_ = c_; c_ = eps[x + 1];
				if (Mask[x] == 1) q[x] -= 0.25 * (c_ - a_) * (c - a) * fjc * fjc / C;
			}
		}
		for (x = fjc; x < MX + fjc; x++)
		if (Mask[x] == 1) {
			g[x] = -q[x];
		}
		return;
	}
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
NAMICS_DBG("remove_bounds in LGrad1 " << std::endl);	int k;
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
NAMICS_DBG("set_bounds in LGrad1 " << std::endl);	int k;
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
NAMICS_DBG("set_bounds in LGrad1 " << std::endl);	int k=0;
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
NAMICS_DBG("set_M_bounds in LGrad1 " << std::endl); //set mirror bounds
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
NAMICS_DBG("remove_bounds in LGrad1 " << std::endl);int k;
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
NAMICS_DBG("set_bounds in LGrad1 " << std::endl);	int k=0;
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

Real LGrad1::ComputeGN(Real* G){
	return WeightedSum(G);
}

void LGrad1::Initiate(Real* G,Real* Gz){
	std::copy_n(Gz, M, G);
}
