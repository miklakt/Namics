#ifndef TOOLS_H
#define TOOLS_H
#include "namics.h"
#include <cassert>
#include <initializer_list>
#include <span>
#include <algorithm>
#include <utility>

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

template<typename F>
inline void for_each_box_site(std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY, F&& apply) {
	assert(Bx.size() >= static_cast<size_t>(n_box));
	assert(By.size() >= static_cast<size_t>(n_box));
	assert(Bz.size() >= static_cast<size_t>(n_box));
	auto project = [](int coord, int limit, int stride) {
		return (coord > limit) ? (coord - limit) * stride : coord * stride;
	};
	for (int p=0; p<n_box; ++p) {
		const int base = p * M;
		const int Bxp = Bx[p];
		const int Byp = By[p];
		const int Bzp = Bz[p];
		for (int i=1; i<=Mx; ++i) {
			const int local_x = i * jx;
			const int pos_x = project(Bxp + i, MX, JX);
			for (int j=1; j<=My; ++j) {
				const int local_y = j * jy;
				const int pos_y = project(Byp + j, MY, JY);
				for (int k=1; k<=Mz; ++k) {
					const int local = base + local_x + local_y + k;
					const int pos_z = project(Bzp + k, MZ, 1);
					apply(p, local, pos_x + pos_y + pos_z);
				}
			}
		}
	}
}

namespace tools {

template<typename T>
void DistributeG1(std::span<const T> G1, std::span<Real> g1, std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	assert(Bx.size() >= static_cast<size_t>(n_box));
	assert(By.size() >= static_cast<size_t>(n_box));
	assert(Bz.size() >= static_cast<size_t>(n_box));
	assert(g1.size() >= static_cast<size_t>(n_box * M));
	for_each_box_site(Bx, By, Bz, M, n_box, Mx, My, Mz, MX, MY, MZ, jx, jy, JX, JY,
	                  [&g1, &G1](int, int local, int global) {
		g1[local] = G1[global];
	});
}

template<typename T>
void CollectPhi(std::span<T> phi, std::span<const Real> GN, std::span<const Real> rho, std::span<const int> Bx, std::span<const int> By, std::span<const int> Bz, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	assert(Bx.size() >= static_cast<size_t>(n_box));
	assert(By.size() >= static_cast<size_t>(n_box));
	assert(Bz.size() >= static_cast<size_t>(n_box));
	assert(GN.size() >= static_cast<size_t>(n_box));
	assert(rho.size() >= static_cast<size_t>(n_box * M));
	for_each_box_site(Bx, By, Bz, M, n_box, Mx, My, Mz, MX, MY, MZ, jx, jy, JX, JY,
	                  [&phi, &GN, &rho](int p, int local, int global) {
		phi[global] += rho[local] / GN[p];
	});
}

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

template <typename T>
inline void add_terms(T* dst, int n, std::initializer_list<std::pair<const T*, T>> terms) {
	for (const auto& [src, scale] : terms) add_shifted(dst, src, n, scale);
}

template <typename T>
inline void add_from_source(const T* src, int n, std::initializer_list<std::pair<T*, T>> targets) {
	for (const auto& [dst, scale] : targets) add_shifted(dst, src, n, scale);
}

#endif
