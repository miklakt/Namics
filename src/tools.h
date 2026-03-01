#ifndef TOOLSxH
#define TOOLSxH
#include "namics.h"
#include <numeric>

extern Real* SUM_RESULT;



#include "tools_host.h"


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
  Real Sum=0;
	for (int i=0; i<M; i++) Sum+=H[i];
  return Sum;
}

template <typename T>
inline Real H_Dot(T* A, T* B, int M) {
  Real Sum=0;
	for (int i = 0; i<M; i++)
    Sum += A[i]*B[i];
  return Sum;
}

template <typename T>
inline void H_Invert(T* KSAM, T* MASK, int M) {
  std::transform(MASK, MASK + M, KSAM, KSAM, [](Real A, Real B) {if (A==0) return 1.0; else return 0.0;});
}

template<typename T>
void H_PutAlpha(T *g, T *phitot, T *phi_side, T chi, T phibulk, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}

Real pythag(Real, Real);
int svdcmp(Real**, int, int, Real*, Real**);

#endif
