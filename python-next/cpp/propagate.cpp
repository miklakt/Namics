#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <Python.h>

// Python owns every buffer. Lattice neighbors and weights are [site, width],
// segment weights and densities are [type, site], and forward is [contour, site].
extern "C" int namics_linear_density(const double* weights, const std::int32_t* chain_type,
    const std::int32_t* neighbors, const double* transitions, int width,
    int types, int sites, int length, double bulk, double* forward,
    double* segment_density, double* molecule_density, double* scratch) {
    if (!weights || !chain_type || !neighbors || !transitions || !forward ||
        !segment_density || !molecule_density || !scratch || width < 1 ||
        types < 1 || sites < 1 || length < 1) return 1;
    for (int s = 0; s < length; ++s)
        if (chain_type[s] < 0 || chain_type[s] >= types) return 2;
    for (std::size_t i = 0; i < static_cast<std::size_t>(sites) * width; ++i)
        if (neighbors[i] < 0 || neighbors[i] >= sites) return 3;
    auto step = [sites, width, neighbors, transitions](const double* in, double* out) {
        for (int i = 0; i < sites; ++i) {
            double value = 0;
            for (int j = 0; j < width; ++j) {
                const auto k = static_cast<std::size_t>(i) * width + j;
                value += transitions[k] * in[neighbors[k]];
            }
            out[i] = value;
        }
    };
    std::copy_n(weights + static_cast<std::size_t>(chain_type[0]) * sites, sites, forward);
    for (int s = 1; s < length; ++s) {
        double* next = forward + static_cast<std::size_t>(s) * sites;
        const double* weight = weights + static_cast<std::size_t>(chain_type[s]) * sites;
        step(next - sites, next);
        for (int i = 0; i < sites; ++i) next[i] *= weight[i];
    }
    double* back = scratch;
    double* next = scratch + sites;
    std::fill_n(back, sites, 1.0);
    for (int s = length - 1; s >= 0; --s) {
        const double* f = forward + static_cast<std::size_t>(s) * sites;
        const double* weight = weights + static_cast<std::size_t>(chain_type[s]) * sites;
        double* out = segment_density + static_cast<std::size_t>(chain_type[s]) * sites;
        for (int i = 0; i < sites; ++i) {
            const double value = (bulk / length) * f[i] * back[i];
            out[i] += value;
            molecule_density[i] += value;
            back[i] *= weight[i];
        }
        if (s) { step(back, next); std::swap(back, next); }
    }
    return 0;
}
static PyModuleDef module = {PyModuleDef_HEAD_INIT, "_namics_kernels", nullptr, -1, nullptr, nullptr, nullptr, nullptr, nullptr};
PyMODINIT_FUNC PyInit__namics_kernels() { return PyModule_Create(&module); }
