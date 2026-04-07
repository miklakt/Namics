#include <iostream>
#include <string>
#include "LGrad3.h"
#include "tools.h"

LGrad3::LGrad3(const Input& In_,const std::string& name_): Lattice(In_,name_) {}

LGrad3::~LGrad3() {
NAMICS_DBG("LGrad3 destructor " << std::endl);}

void LGrad3::ComputeLambdas() {
}

bool LGrad3::CheckLatticeInput(const ParameterStore& parameters) {
	bool success = RejectScalarBoundsInMultiD(parameters);
	success = ReadScaledDimension(parameters, "n_layers_x", MX, 1, "In 'lat' the parameter 'n_layers_x' is required", "n_layers_x out of bounds, currently: 1..1e6; Problem terminated") && success;
	success = ReadScaledDimension(parameters, "n_layers_y", MY, 1, "In 'lat' the parameter 'n_layers_y' is required", "n_layers_y out of bounds, currently: 1..1e6; Problem terminated") && success;
	success = ReadScaledDimension(parameters, "n_layers_z", MZ, 1, "In 'lat' the parameter 'n_layers_z' is required", "n_layers_z out of bounds, currently: 1..1e6; Problem terminated") && success;

	success = ReadBoundaryCondition(parameters, "lowerbound_x", 0, {"mirror", "periodic"}, "for 'lowerbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ") && success;
	success = ReadBoundaryCondition(parameters, "upperbound_x", 3, {"mirror", "periodic"}, "for 'upperbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ") && success;
	success = ReadBoundaryCondition(parameters, "lowerbound_y", 1, {"mirror", "periodic"}, "for 'lowerbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ") && success;
	success = ReadBoundaryCondition(parameters, "upperbound_y", 4, {"mirror", "periodic"}, "for 'upperbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ") && success;
	success = ReadBoundaryCondition(parameters, "lowerbound_z", 2, {"mirror", "periodic"}, "for 'lowerbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ") && success;
	success = ReadBoundaryCondition(parameters, "upperbound_z", 5, {"mirror", "periodic"}, "for 'upperbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ") && success;
	success = CheckPeriodicPair(1, 4, ("In y-direction the boundary conditions do not match:" + BC[1] + " and " + BC[4]).c_str()) && success;
	success = CheckPeriodicPair(2, 5, ("In z-direction the boundary conditions do not match:" + BC[2] + " and " + BC[5]).c_str()) && success;
	return success;
}

bool LGrad3::PutM() {
NAMICS_DBG("PutM in LGrad3 " << std::endl);	bool success=true;
	volume = MX*MY*MZ;
	JX=(MZ+2*fjc)*(MY+2*fjc); JY=MZ+2*fjc; JZ=1; M = (MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc);

	Accesible_volume=volume;
	return success;
}

Real LGrad3:: Moment(Real* X,Real Xb, int n) {
	(void)n;
	(void)Xb;
	(void)X;
NAMICS_DBG("Moment in LGrad3 " << std::endl);	Real Result=0;
	return Result/fjc;
}

Real LGrad3::WeightedSum(Real* X){
NAMICS_DBG("weighted sum in LGrad3 " << std::endl);
	remove_bounds(X);
	return std::accumulate(X, X + M, Real{0});
}

void LGrad3::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
NAMICS_DBG(" Side in LGrad3 " << std::endl);	if (ignore_sites) {
		std::copy_n(X, M, X_side); return;
	}
	std::fill_n(X_side, M, 0);//set_bounds(X);
	Real Two=2.0;
	Real Four=4.0;
	if (!stencil_full) {
		add_shifted(X_side + JX, X, M - JX);
		add_shifted(X_side, X + JX, M - JX);
		add_shifted(X_side + JY, X, M - JY);
		add_shifted(X_side, X + JY, M - JY);
		add_shifted(X_side + 1, X, M - 1);
		add_shifted(X_side, X + 1, M - 1);
		scale_span(X_side, M, lattice_type == simple_cubic ? Four : Two);

		add_shifted(X_side + JX + JY, X, M - JX - JY);
		add_shifted(X_side, X + JX + JY, M - JX - JY);
		add_shifted(X_side + JY, X + JX, M - JY - JX);
		add_shifted(X_side + JX, X + JY, M - JY - JX);
		add_shifted(X_side + JX + 1, X, M - JX - 1);
		add_shifted(X_side, X + JX + 1, M - JX - 1);
		add_shifted(X_side + JX, X + 1, M - JX);
		add_shifted(X_side + 1, X + JX, M - JX);
		add_shifted(X_side + JY + 1, X, M - JY - 1);
		add_shifted(X_side, X + JY + 1, M - JX - 1);
		add_shifted(X_side + JY, X + 1, M - JY);
		add_shifted(X_side + 1, X + JY, M - JY);
		scale_span(X_side, M, lattice_type == simple_cubic ? Four : Two);

		add_shifted(X_side + JX + JY + 1, X, M - JX - JY - 1);
		add_shifted(X_side, X + JX + JY + 1, M - JX - JY - 1);
		add_shifted(X_side + JX + JY, X + 1, M - JX - JY - 1);
		add_shifted(X_side + 1, X + JX + JY, M - JX - JY - 1);
		add_shifted(X_side + JX + 1, X + JY, M - JX - JY - 1);
		add_shifted(X_side + JY, X + JX + 1, M - JX - JY - 1);
		add_shifted(X_side + JY + 1, X + JX, M - JX - JY - 1);
		add_shifted(X_side + JX, X + JY + 1, M - JX - JY - 1);
		scale_span(X_side, M, lattice_type == simple_cubic ? Real(1.0 / 152.0) : Real(1.0 / 56.0));
	} else {
		if (lattice_type==simple_cubic) {
			Real C=1.0/6.0;
			for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
			for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
			for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (X_side+JZ)[__i] += (X)[__i];
			for (int __i = 0; __i < (M-JZ); ++__i) (X_side)[__i] += (X+JZ)[__i];
	 		scale_span(X_side, M, C);
		} else { //hexagonal
			if (fjc==1) {
				Real Two=2.0;
				Real C=1.0/40.0;
				add_shifted(X_side, X, M, Two);
				for (int __i = 0; __i < (M-JX); ++__i) (X_side+JX)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (X_side)[__i] += (X+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (X_side+JY)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (X_side)[__i] += (X+JY)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side+1)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (X_side)[__i] += (X+1)[__i];
				scale_span(X_side, M, Two);

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
				scale_span(X_side, M, Two);

				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+JY+1)[__i] += (X)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side)[__i] += (X+JX+JY+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+JY)[__i] += (X+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+1)[__i] += (X+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX+1)[__i] += (X+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JY)[__i] += (X+JX+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JY+1)[__i] += (X+JX)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (X_side+JX)[__i] += (X+JY+1)[__i];

				scale_span(X_side, M, C);


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
					if (block !=3) scale_span(X_side, M, Two); else scale_span(X_side, M, C);

				}
			}
		}
	}
}

void LGrad3::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
NAMICS_DBG(" propagate in LGrad3 " << std::endl); Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
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
		if (lattice_type == simple_cubic) scale_span(gs, M, 4.0);
		else scale_span(gs, M, 2.0);
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
		if (lattice_type == simple_cubic) scale_span(gs, M, 4.0);
		else scale_span(gs, M, 2.0);
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_+JY_+1)[__i] += (gs_1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs)[__i] += (gs_1+JX_+JY_+1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_+JY_)[__i] += (gs_1+1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+1)[__i] += (gs_1+JX_+JY_)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_+1)[__i] += (gs_1+JY_)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JY_)[__i] += (gs_1+JX_+1)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JY_+1)[__i] += (gs_1+JX_)[__i];
		for (int __i = 0; __i < (M-JX_-JY_-1); ++__i) (gs+JX_)[__i] += (gs_1+JY_+1)[__i];
		if (lattice_type == simple_cubic) scale_span(gs, M, 1.0 / 152.0);
		else scale_span(gs, M, 1.0 / 56.0);
		std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });
	} else {
		if (lattice_type==simple_cubic) {
			for (int __i = 0; __i < (M-JX_); ++__i) (gs+JX_)[__i] += (gs_1)[__i];
			for (int __i = 0; __i < (M-JX_); ++__i) (gs)[__i] += (gs_1+JX_)[__i];
			for (int __i = 0; __i < (M-JY_); ++__i) (gs+JY_)[__i] += (gs_1)[__i];
			for (int __i = 0; __i < (M-JY_); ++__i) (gs)[__i] += (gs_1+JY_)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (gs_1)[__i];
			for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (gs_1+1)[__i];
			scale_span(gs, M, 1.0 / 6.0);
			std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });
		} else { //hexagonal
			if (fjc==1) {
				Real Two=2.0;
				Real C=1.0/40.0;
				add_shifted(gs, gs_1, M, Two);
				for (int __i = 0; __i < (M-JX); ++__i) (gs+JX)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX); ++__i) (gs)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs+JY)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JY); ++__i) (gs)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (gs+1)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-1); ++__i) (gs)[__i] += (gs_1+1)[__i];
				scale_span(gs, M, Two);

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
				scale_span(gs, M, Two);

				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX+JY+1)[__i] += (gs_1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs)[__i] += (gs_1+JX+JY+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX+JY)[__i] += (gs_1+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+1)[__i] += (gs_1+JX+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX+1)[__i] += (gs_1+JY)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JY)[__i] += (gs_1+JX+1)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JY+1)[__i] += (gs_1+JX)[__i];
				for (int __i = 0; __i < (M-JX-JY-1); ++__i) (gs+JX)[__i] += (gs_1+JY+1)[__i];

				scale_span(gs, M, C);
				std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });


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
					if (block !=3) scale_span(gs, M, 2.0); else scale_span(gs, M, 1.0 / 8.0 / ((FJC-2)*(FJC-2)*(FJC-2)+3*(FJC-2)*(FJC-2)+3*(FJC-2)+1));

				}
				std::transform(gs, gs + M, G1, gs, [](auto a, auto b) { return a * b; });


			}
		}
	}
}

void LGrad3::UpdateEE(Real* EE, Real* psi, Real* E) {
	(void)E;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);

	std::fill_n(EE, M, 0);
	for (int __i = 0; __i < (M-2); ++__i) (EE+1)[__i] += std::pow((psi)[__i]-(psi+1)[__i],2) + std::pow((psi+1)[__i]-(psi+2)[__i],2);
	for (int __i = 0; __i < (M-2*JX); ++__i) (EE+JX)[__i] += std::pow((psi)[__i]-(psi+JX)[__i],2) + std::pow((psi+JX)[__i]-(psi+2*JX)[__i],2);
	for (int __i = 0; __i < (M-2*JY); ++__i) (EE+JY)[__i] += std::pow((psi)[__i]-(psi+JY)[__i],2) + std::pow((psi+JY)[__i]-(psi+2*JY)[__i],2);
	scale_span(EE, M, pf);


}


void LGrad3::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, Real* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	(void)grad_epsilon;
	int x, y, z;

	Real epsZplus, epsZmin, epsXplus, epsXmin, epsYplus, epsYmin;
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
NAMICS_DBG("set_bounds_x (shift in y,z ) in LGrad3 " << std::endl);	int y,z;
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
NAMICS_DBG("set_bounds_y (shift in x,z) in LGrad3 " << std::endl);	int x,z;
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
NAMICS_DBG("set_bounds_z shift (x,y) in LGrad3 " << std::endl);	int x,y;
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
NAMICS_DBG("set_bounds_x (shift yz) in LGrad3 " << std::endl);	int y,z;
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
NAMICS_DBG("set_bounds_y (shift x,z) in LGrad3 " << std::endl);	int x,z;
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
NAMICS_DBG("set_bounds_z (shift xy) in LGrad3 " << std::endl);	int x,y;
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
NAMICS_DBG("remove_bounds in LGrad3 " << std::endl);	int x,y,z;
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
NAMICS_DBG("set_bounds in LGrad3 " << std::endl);	int x,y,z;
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
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(std::span<Real>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {
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
				for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[x*JX+(2*fjc-1-k)*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+(MY+fjc-k-1)*JY+z*JZ];
			}
		}
	}
}

void LGrad3::remove_bounds(int *X){
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
				for (z=0; z<MZ+2*fjc; z++) for (x=0; x<MX+2*fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = 0;
			}
		}
	}
}

void LGrad3::set_bounds(int* X){
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(std::span<int>(X+i*m[k], static_cast<size_t>(m[k])),jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {

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

Real LGrad3::ComputeGN(Real* G, int M){
	return WeightedSum(G);
}

void LGrad3::Initiate(Real* G,Real* Gz){
	std::copy_n(Gz, M, G);
}

bool LGrad3:: PutMask(Real* MASK,std::vector<int>px,std::vector<int>py,std::vector<int>pz,int R){
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
	}

	return success;
}
