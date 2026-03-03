#ifndef HOST_ONLY_TOOLSxH
#define HOST_ONLY_TOOLSxH
#include "namics.h"
#include <cassert>
#include <span>
#include <algorithm>
#include <numeric>

template <typename T>
inline void H_Cp(Real* P, T* A, int M) {
	std::copy(A, A + M, P);
}

template <typename T>
inline void H_Zero(T* H_P, int M) {
	std::fill(H_P, H_P + M, 0);
}

template <typename T>
inline Real H_Sum(T* H, int M) {
	Real Sum = 0;
	for (int i = 0; i < M; i++) Sum += H[i];
	return Sum;
}

template <typename T>
inline Real H_Dot(T* A, T* B, int M) {
	Real Sum = 0;
	for (int i = 0; i < M; i++) Sum += A[i] * B[i];
	return Sum;
}

template <typename T>
inline void H_Invert(T* KSAM, T* MASK, int M) {
	std::transform(MASK, MASK + M, KSAM, KSAM, [](Real A, Real B) { if (A == 0) return 1.0; else return 0.0; });
}

template<typename T>
void H_PutAlpha(T *g, T *phitot, T *phi_side, T chi, T phibulk, int M) {
	for (int i = 0; i < M; i++) if (phitot[i] > 0) g[i] = g[i] - chi * (phi_side[i] / phitot[i] - phibulk);
}

template <typename T>
void bx(std::span<T> P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy, int by1=0, int bz1=0, bool corner=false)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	int jy_by1=jy*by1;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
		if (corner) {
			if (y==0 && z==0) {
				P[i]=P[bx1_jx+i+jy_by1+bz1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1+bz1];
			}
			if (y==My-1 && z==0) {
				P[i]=P[bx1_jx+i-jy_by1+bz1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1+bz1];
			}
			if (y==0 && z==Mz-1) {
				P[i]=P[bx1_jx+i+jy_by1-bz1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1-bz1];
			}
			if (y==My-1 && z==Mz-1) {
				P[i]=P[bx1_jx+i-jy_by1-bz1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1-bz1];
			}
			if (y==0) {
				P[i]=P[bx1_jx+i+jy_by1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1];
			}
			if (y==My-1) {
				P[i]=P[bx1_jx+i-jy_by1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1];
			}
			if (z==0) {
				P[i]=P[bx1_jx+i+bz1];
				P[jx_mmx+i]=P[jx_bxm+i+bz1];
			}
			if (z==Mz-1) {
				P[i]=P[bx1_jx+i-bz1];
				P[jx_mmx+i]=P[jx_bxm+i-bz1];
			}
		}
	}
}

template<typename T>
void b_x(std::span<T> P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	(void)bxm;
	(void)bx1;
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}

template<typename T>
void by(std::span<T> P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy, int bz1=0, bool corner=false)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
		if (corner) {
			if (z==0) {
				P[i]=P[jy_by1+i+bz1];
				P[jy_mmy+i]=P[jy_bym+i+bz1];
			}
			if (z==Mz-1) {
				P[i]=P[jy_by1+i-bz1];
				P[jy_mmy+i]=P[jy_bym+i-bz1];
			}
		}
	}
}

template<typename T>
void b_y(std::span<T> P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	(void)bym;
	(void)by1;
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}

template<typename T>
void bz(std::span<T> P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy, int bx1=0, bool corner=false)   {
	int i;
	int bx1_jx=bx1*jx;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
		if (corner) {
			if (x==0) {
				P[i]=P[i+bz1+bx1_jx];
				P[i+mmz]=P[i+bzm+bx1_jx];
			}
			if (x==Mx-1) {
				P[i]=P[i+bz1-bx1_jx];
				P[i+mmz]=P[i+bzm-bx1_jx];
			}
		}
	}
}

template<typename T>
void b_z(std::span<T> P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	(void)bzm;
	(void)bz1;
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}

template <typename T>
inline void SetBoundaries(std::span<T> P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz, bool corners=false) {
  bx(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy, by1, bz1, corners);
  by(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy, bz1, corners);
  bz(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy, bx1, corners);
}

template <typename T>
inline void RemoveBoundaries(std::span<T> P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  b_x(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy);
  b_y(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy);
  b_z(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy);
}

template<typename T>
void Xr_times_ci(int posi, int k_diis, int k, int m, std::span<T> x, std::span<const T> xR, std::span<const T> Ci) {
	const int nvar = static_cast<int>(x.size());
	assert(xR.size() >= static_cast<size_t>(m * nvar));
	assert(Ci.size() >= static_cast<size_t>(k_diis));
	for (int __i = 0; __i < nvar; ++__i) {
		x[__i] += Ci[0] * xR[posi*nvar + __i];
	}

	for (int i=1; i<k_diis; i++) {
		posi = k-k_diis+1+i;
	    	if (posi<0) {
	      		posi +=m;
		}
		for (int __j = 0; __j < nvar; ++__j) {
			x[__j] += Ci[i] * xR[posi*nvar + __j];
		}
	}
}

template<typename T>
inline void Xr_times_ci(int posi, int k_diis, int k, int m, int nvar, T* x, const T* xR, const T* Ci) {
	assert(nvar >= 0);
	Xr_times_ci(posi, k_diis, k, m,
	            std::span<T>(x, static_cast<size_t>(nvar)),
	            std::span<const T>(xR, static_cast<size_t>(m * nvar)),
	            std::span<const T>(Ci, static_cast<size_t>(k_diis)));
}

namespace tools {

template<typename T>
void DistributeG1(std::span<const T> G1, std::span<Real> g1, std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	(void)MM;
	assert(Bx.size() >= static_cast<size_t>(n_box));
	assert(By.size() >= static_cast<size_t>(n_box));
	assert(Bz.size() >= static_cast<size_t>(n_box));
	assert(g1.size() >= static_cast<size_t>(n_box * M));
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p];
		for (int i=1; i<Mx+1; i++) { ii+=jx; jj=0; if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					g1[pos_l+ii+jj+kk]=G1[pos_x+pos_y+pos_z];
				}
			}
		}
	}
}

template<typename T>
void CollectPhi(std::span<T> phi, std::span<const Real> GN, std::span<const Real> rho, std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	(void)MM;
	assert(Bx.size() >= static_cast<size_t>(n_box));
	assert(By.size() >= static_cast<size_t>(n_box));
	assert(Bz.size() >= static_cast<size_t>(n_box));
	assert(GN.size() >= static_cast<size_t>(n_box));
	assert(rho.size() >= static_cast<size_t>(n_box * M));
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	Real Inv_H_GNp;
	int ii=0,jj=0,kk=0;
	for (int p=0; p<n_box; p++) {pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p]; Inv_H_GNp=1.0/GN[p];
		for (int i=1; i<Mx+1; i++) {ii+=jx; jj=0;  if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					phi[pos_x+pos_y+pos_z]+=rho[pos_l+ii+jj+kk]*Inv_H_GNp;
				}
			}
		}
	}
}

}

template<typename T, typename D>
void OverwriteC(std::span<T> P, std::span<const D> Mask, T C) {
	assert(P.size() == Mask.size());
	std::transform(P.begin(), P.end(), Mask.begin(), P.begin(),
	               [C](T p, D mask_value) { return (mask_value == 1) ? C : p; });
}

template<typename T, typename D>
inline void OverwriteC(T* P, const D* Mask, T C, int M) {
	assert(M >= 0);
	OverwriteC(std::span<T>(P, static_cast<size_t>(M)),
	           std::span<const D>(Mask, static_cast<size_t>(M)), C);
}

template<typename T, typename D>
void OverwriteA(std::span<T> P, std::span<const D> Mask, std::span<const T> A) {
	assert(P.size() == Mask.size());
	assert(P.size() == A.size());
	std::transform(A.begin(), A.end(), Mask.begin(), P.begin(),
	               [](T a, D mask_value) { return (mask_value == 1) ? a : T{0}; });
}

template<typename T, typename D>
inline void OverwriteA(T* P, const D* Mask, const T* A, int M) {
	assert(M >= 0);
	OverwriteA(std::span<T>(P, static_cast<size_t>(M)),
	           std::span<const D>(Mask, static_cast<size_t>(M)),
	           std::span<const T>(A, static_cast<size_t>(M)));
}

#endif
