#ifndef HOST_ONLY_TOOLSxH
#define HOST_ONLY_TOOLSxH

#include <numeric>
#include <algorithm>
#include <functional>
#if __SSE__
	#include <smmintrin.h>
#endif
#include <cmath>


//typedef double Real;
struct saxpy_functor
{
	const Real a;

	saxpy_functor(Real _a) : a(_a) {}

	Real operator()(const Real &x, const Real &y) const
	{
		return a * x + y;
	}
};

struct reverse_minus_functor
{
	reverse_minus_functor() {}

	Real operator()(const Real &x, const Real &y) const
	{
		return y - x;
	}
};

struct norm_functor
{
	const Real a;

	norm_functor(Real _a) : a(_a) {}

	Real operator()(const Real &x) const
	{
		return a * x;
	}
};

struct binary_norm_functor
{
	const Real a;

	binary_norm_functor(Real _a) : a(_a) {}

	Real operator()(const Real &x, const Real &y) const
	{
		return a * x * y;
	}
};

struct order_param_functor
{

	order_param_functor() {}

	Real operator()(const Real &x, const Real &y) const
	{
		return pow(x - y, 2);
	}
};

struct is_negative_functor
{

	const Real tolerance{0};

	is_negative_functor(Real _tolerance = 0) : tolerance(_tolerance) {}

	bool operator()(const Real &x) const
	{
		return x < 0 - tolerance || x > 1 + tolerance;
	}
};

struct is_not_unity_functor
{
	const Real tolerance{0};

	is_not_unity_functor(Real _tolerance = 1e-4) : tolerance(_tolerance) {}

	bool operator()(const Real &x) const
	{
		bool result{0};

		if (x > (1 + tolerance) || x < (1 - tolerance))
			result = 1;

		return result;
	}
};

//typedef long double Real;


template <typename T>
void bx(T* P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy, int by1=0, int bz1=0, bool corner=false)   {
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
void b_x(T *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
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
void by(T *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy, int bz1=0, bool corner=false)   {
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
void b_y(T *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
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
void bz(T *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy, int bx1=0, bool corner=false)   {
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
void b_z(T *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
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
inline void SetBoundaries(T* P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz, bool corners=false) {
  bx(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy, by1, bz1, corners);
  by(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy, bz1, corners);
  bz(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy, bx1, corners);
}

template <typename T>
inline void RemoveBoundaries(T* P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  b_x(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy);
  b_y(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy);
  b_z(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy);
}

//template<typename T>
//	T z = 0.0;
//	T ftmp[2] = { zero, zero };
//	__m128d mres;
//
//			_mm_loadu_pd(&y[2*i])));
//
//		_mm_store_pd(ftmp, mres);
//
//}
//
//			result += x[i] * y[i];
//	}
//}


template<typename T>
void Xr_times_ci(int posi, int k_diis, int k, int m, int nvar, T* x, T* xR, T* Ci) {
	for (int __i = 0; __i < nvar; ++__i) x[__i] += Ci[0] * xR[posi*nvar + __i];

	for (int i=1; i<k_diis; i++) {
		posi = k-k_diis+1+i;
	    	if (posi<0) {
	      		posi +=m;
		}
		for (int __j = 0; __j < nvar; ++__j) x[__j] += Ci[i] * xR[posi*nvar + __j];
	}
}

namespace tools {

template<typename T>
void DistributeG1(T* G1, Real* g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	(void)MM;
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
void CollectPhi(T* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	(void)MM;
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
void OverwriteC(T *P, D *Mask, T C, int M) {
	for (int i=0; i<M; i++) if (Mask[i]==1) P[i]=C; //else P[i]=0;
}

template<typename T, typename D>
void OverwriteA(T *P, D *Mask,T* A,int M) {
	for (int i=0; i<M; i++) if (Mask[i]==1) P[i]=A[i]; else P[i]=0;
}

Real pythag(Real, Real);
int svdcmp(Real**, int, int, Real*, Real**);

#endif
