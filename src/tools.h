#ifndef TOOLS_H
#define TOOLS_H
#include "namics.h"
#include <span>

template<typename F>
inline void for_each_x_face(int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy, int by1, int bz1, bool corners, F&& apply) {
	const int jx_mmx = jx * mmx;
	const int jx_bxm = jx * bxm;
	const int bx1_jx = bx1 * jx;
	const int jy_by1 = jy * by1;
	for (int y=0; y<My; ++y) {
		const int y_shift = corners ? ((y == 0 ? jy_by1 : 0) + (y == My - 1 ? -jy_by1 : 0)) : 0;
		for (int z=0; z<Mz; ++z) {
			const int z_shift = corners ? ((z == 0 ? bz1 : 0) + (z == Mz - 1 ? -bz1 : 0)) : 0;
			const int i = jy * y + z;
			apply(i, jx_mmx + i, bx1_jx + i + y_shift + z_shift, jx_bxm + i + y_shift + z_shift);
		}
	}
}

template<typename F>
inline void for_each_y_face(int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy, int bz1, bool corners, F&& apply) {
	const int jy_mmy = jy * mmy;
	const int jy_bym = jy * bym;
	const int jy_by1 = jy * by1;
	for (int x=0; x<Mx; ++x) {
		for (int z=0; z<Mz; ++z) {
			const int z_shift = corners ? ((z == 0 ? bz1 : 0) + (z == Mz - 1 ? -bz1 : 0)) : 0;
			const int i = jx * x + z;
			apply(i, jy_mmy + i, jy_by1 + i + z_shift, jy_bym + i + z_shift);
		}
	}
}

template<typename F>
inline void for_each_z_face(int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy, int bx1, bool corners, F&& apply) {
	const int bx1_jx = bx1 * jx;
	for (int x=0; x<Mx; ++x) {
		const int x_shift = corners ? ((x == 0 ? bx1_jx : 0) + (x == Mx - 1 ? -bx1_jx : 0)) : 0;
		for (int y=0; y<My; ++y) {
			const int i = jx * x + jy * y;
			apply(i, i + mmz, i + bz1 + x_shift, i + bzm + x_shift);
		}
	}
}

template <typename T>
inline void SetBoundaries(std::span<T> P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz, bool corners=false) {
  for_each_x_face(Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy, by1, bz1, corners,
                  [&P](int dst_low, int dst_high, int src_low, int src_high) {
		P[dst_low] = P[src_low];
		P[dst_high] = P[src_high];
	});
  for_each_y_face(Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy, bz1, corners,
                  [&P](int dst_low, int dst_high, int src_low, int src_high) {
		P[dst_low] = P[src_low];
		P[dst_high] = P[src_high];
	});
  for_each_z_face(Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy, bx1, corners,
                  [&P](int dst_low, int dst_high, int src_low, int src_high) {
		P[dst_low] = P[src_low];
		P[dst_high] = P[src_high];
	});
}

template <typename T>
inline void RemoveBoundaries(std::span<T> P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  for_each_x_face(Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy, 0, 0, false,
                  [&P](int dst_low, int dst_high, int, int) {
		P[dst_low] = T{};
		P[dst_high] = T{};
	});
  for_each_y_face(Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy, 0, false,
                  [&P](int dst_low, int dst_high, int, int) {
		P[dst_low] = T{};
		P[dst_high] = T{};
	});
  for_each_z_face(Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy, 0, false,
                  [&P](int dst_low, int dst_high, int, int) {
		P[dst_low] = T{};
		P[dst_high] = T{};
	});
}

template <typename T>
inline void add_shifted(T* dst, const T* src, int n, T scale = T{1}) {
	for (int i = 0; i < n; ++i) dst[i] += scale * src[i];
}

template <typename T>
inline void add_weighted(T* dst, const T* src, const T* weights, int n, T scale = T{1}) {
	for (int i = 0; i < n; ++i) dst[i] += scale * src[i] * weights[i];
}

template <typename T>
inline void scale_span(T* dst, int n, T scale) {
	for (int i = 0; i < n; ++i) dst[i] *= scale;
}

#endif
