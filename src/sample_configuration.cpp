#include "sample_configuration.h"

#include <random>

SampleConfiguration::SampleConfiguration(std::string name_, const Molecule* molecule_, uint64_t seed_, int size_)
	: name(std::move(name_)), molecule_name(molecule_->name), seed(seed_), size(size_), molecule(molecule_) {}

bool SampleConfiguration::Generate() {
	const Lattice& lat = *molecule->lat;
	const int M = lat.M;
	// One engine serves the whole batch: a fixed seed reproduces both every sample and their ordering.
	std::mt19937_64 random(seed);
	std::uniform_real_distribution<Real> uniform(0, 1);
	std::normal_distribution<Real> normal;
	ParameterStore configurations = ParameterStore::array();
	// Partition Euclidean radius so every curved shell has exactly the SCFT measure L.
	// Its volume-median radius is a deterministic representative of that sampled layer.
	std::vector<Real> layer_coordinates(lat.M);
	Real edge = lat.offset_first_layer / lat.fjc;
	for (int layer = lat.fjc; layer < lat.MX + lat.fjc; ++layer) {
		Real outer;
		if (lat.geometry == "cylindrical") {
			outer = std::sqrt(edge * edge + lat.L[layer] / PIE);
			layer_coordinates[layer] = std::sqrt((edge * edge + outer * outer) / 2);
		} else if (lat.geometry == "spherical") {
			outer = std::cbrt(edge * edge * edge + 3 * lat.L[layer] / (4 * PIE));
			layer_coordinates[layer] = std::cbrt((edge * edge * edge + outer * outer * outer) / 2);
		} else {
			outer = edge + lat.L[layer];
			layer_coordinates[layer] = (edge + outer) / 2;
		}
		edge = outer;
	}
	const auto layer_coordinate = [&](int layer) { return layer_coordinates[layer]; };

	for (int sample = 0; sample < size; ++sample) {
		std::vector<int> layers(molecule->chainlength);
		std::vector<std::array<Real, 3>> coordinates(molecule->chainlength);
		std::vector<Real> weights;
		auto root = std::span<const Real>(molecule->q_forward).first(M);
		// q_forward already sums every continuation below the root. L is required because one radial layer
		// represents L[layer] equivalent sites; omitting it would over-sample the small inner shells.
		for (int layer = lat.fjc; layer < lat.MX + lat.fjc; ++layer) {
			weights.push_back(lat.L[layer] * root[layer]);
		}
		if (std::any_of(weights.begin(), weights.end(), [](Real weight) { return weight < 0 || !std::isfinite(weight); }) ||
		    std::none_of(weights.begin(), weights.end(), [](Real weight) { return weight > 0; })) {
			std::cout << "Unable to sample '" << name << "': the root propagator has no positive finite weight." << std::endl;
			return false;
		}
		layers[0] = lat.fjc + std::discrete_distribution<int>(weights.begin(), weights.end())(random);
		const Real root_position = layer_coordinate(layers[0]);
		// A one-gradient SCF solution fixes only x, cylindrical radius, or spherical radius. The remaining
		// coordinates are symmetry-equivalent, so choose an unbiased orientation (and fix irrelevant translations).
		if (lat.geometry == "planar") {
			coordinates[0] = {root_position, 0, 0};
		} else if (lat.geometry == "cylindrical") {
			const Real angle = 2 * PIE * uniform(random);
			coordinates[0] = {root_position * std::cos(angle), root_position * std::sin(angle), 0};
		} else {
			std::array<Real, 3> direction{normal(random), normal(random), normal(random)};
			const Real length = std::hypot(direction[0], direction[1], direction[2]);
			for (int axis = 0; axis < 3; ++axis) coordinates[0][axis] = root_position * direction[axis] / length;
		}

		for (int parent = 0; parent < molecule->chainlength; ++parent) {
			for (int child : molecule->segment_path[parent].children) {
				std::vector<int> candidate_layers;
				std::vector<std::pair<int, Real>> transitions;
				// These are the same reduced transition coefficients used by LGrad1::propagate. Keeping the
				// coefficients here, rather than inventing Cartesian moves, preserves the solved lattice measure.
				if (lat.fjc == 1) {
					transitions = {{layers[parent] - 1, lat.lambda_1[layers[parent]]},
					               {layers[parent], lat.lambda0[layers[parent]]},
					               {layers[parent] + 1, lat.lambda1[layers[parent]]}};
				} else {
					transitions.emplace_back(layers[parent], lat.LAMBDA[lat.fjc * M + layers[parent]]);
					for (int j = 0; j < lat.FJC / 2; ++j) {
						const int shift = lat.fjc - j;
						transitions.emplace_back(layers[parent] - shift, lat.LAMBDA[j * M + layers[parent]]);
						transitions.emplace_back(layers[parent] + shift, lat.LAMBDA[(lat.FJC - j - 1) * M + layers[parent]]);
					}
				}
				weights.clear();
				auto q_child = std::span<const Real>(molecule->q_forward).subspan(static_cast<size_t>(child) * M, M);
				for (auto [layer, transition] : transitions) {
					// Apply the lattice boundary map to proposed ghost layers. Repeated mapped candidates are
					// intentionally retained: discrete_distribution then adds their separate transition masses.
					if (layer < lat.fjc) layer = lat.fjc == 1 ? lat.BX1 : lat.B_X1[layer];
					else if (layer >= lat.MX + lat.fjc) layer = lat.fjc == 1 ? lat.BXM : lat.B_XM[layer - lat.MX - lat.fjc];
					if (layer < lat.fjc || layer >= lat.MX + lat.fjc) continue;
					candidate_layers.push_back(layer);
					// Exact PBCG conditional: bare parent-to-child transition times the partition sum of
					// every possible continuation rooted at that child layer.
					weights.push_back(transition * q_child[layer]);
				}
				if (std::any_of(weights.begin(), weights.end(), [](Real weight) { return weight < 0 || !std::isfinite(weight); }) ||
				    std::none_of(weights.begin(), weights.end(), [](Real weight) { return weight > 0; })) {
					std::cout << "Unable to sample '" << name << "': a child propagator has no positive finite weight." << std::endl;
					return false;
				}
				layers[child] = candidate_layers[std::discrete_distribution<size_t>(weights.begin(), weights.end())(random)];
				const Real target = layer_coordinate(layers[child]);
				const auto& origin = coordinates[parent];
				auto& point = coordinates[child];

				// This second stage is a conditional lift of an already sampled layer path. It preserves the
				// projected SCF statistics but does not claim to reconstruct angular correlations absent in 1D.
				if (lat.geometry == "planar") {
					Real x = target;
					if (lat.BC[0] == "periodic") {
						// Select the periodic image reachable from the parent without changing the output layer.
						const Real period = static_cast<Real>(lat.MX) / lat.fjc;
						x += std::round((origin[0] - x) / period) * period;
					}
					const Real dx = x - origin[0];
					const Real lateral = std::sqrt(std::max(Real{0}, 1 - dx * dx));
					const Real angle = 2 * PIE * uniform(random);
					point = {x, origin[1] + lateral * std::cos(angle), origin[2] + lateral * std::sin(angle)};
				} else if (lat.geometry == "cylindrical") {
					// Pick an admissible transverse chord between the two sampled radii; the remaining bond
					// length is placed along the invariant cylinder axis, with no preferred sign or azimuth.
					const Real radius = std::hypot(origin[0], origin[1]);
					const Real low = std::abs(radius - target);
					const Real high = std::min(Real{1}, radius + target);
					const Real transverse = high > low ? low + (high - low) * uniform(random) : low;
					Real angle = 2 * PIE * uniform(random);
					if (radius * target > 0) {
						const Real cosine = std::clamp((radius * radius + target * target - transverse * transverse) / (2 * radius * target), Real{-1}, Real{1});
						angle = std::atan2(origin[1], origin[0]) + (uniform(random) < 0.5 ? -1 : 1) * std::acos(cosine);
					}
					point = {target * std::cos(angle), target * std::sin(angle),
					         origin[2] + (uniform(random) < 0.5 ? -1 : 1) * std::sqrt(std::max(Real{0}, 1 - transverse * transverse))};
				} else {
					// The law of cosines fixes the radial component of a unit bond. A random tangent direction
					// supplies the unresolved spherical angle without changing either sampled shell.
					const Real radius = std::hypot(origin[0], origin[1], origin[2]);
					std::array<Real, 3> direction{normal(random), normal(random), normal(random)};
					Real length = std::hypot(direction[0], direction[1], direction[2]);
					for (Real& value : direction) value /= length;
					if (radius * target > 0) {
						std::array<Real, 3> radial{origin[0] / radius, origin[1] / radius, origin[2] / radius};
						Real projection = direction[0] * radial[0] + direction[1] * radial[1] + direction[2] * radial[2];
						for (int axis = 0; axis < 3; ++axis) direction[axis] -= projection * radial[axis];
						length = std::hypot(direction[0], direction[1], direction[2]);
						if (length == 0) {
							direction = std::abs(radial[0]) < 0.9 ? std::array<Real, 3>{0, -radial[2], radial[1]} : std::array<Real, 3>{-radial[1], radial[0], 0};
							length = std::hypot(direction[0], direction[1], direction[2]);
						}
						const Real cosine = std::clamp((radius * radius + target * target - 1) / (2 * radius * target), Real{-1}, Real{1});
						for (int axis = 0; axis < 3; ++axis) direction[axis] = cosine * radial[axis] + std::sqrt(1 - cosine * cosine) * direction[axis] / length;
					}
					for (int axis = 0; axis < 3; ++axis) point[axis] = target * direction[axis];
				}
			}
		}

		ParameterStore chain = ParameterStore::array();
		for (const auto& point : coordinates) chain.push_back({point[0], point[1], point[2]});
		configurations.push_back(std::move(chain));
	}
	xyz = size == 1 ? std::move(configurations[0]) : std::move(configurations);
	return true;
}

namespace sample_configuration {

bool Load(const Input& input, Lattice* lattice, std::span<const std::unique_ptr<Molecule>> molecules, int start, std::vector<SampleConfiguration>& samples) {
	const auto& problem = input.Start(start);
	const auto section = problem.find("sample");
	if (section == problem.end()) return true;
	if (lattice->gradients != 1) {
		std::cout << "Sample configurations currently require a one-gradient lattice." << std::endl;
		return false;
	}
	for (auto request = section->begin(); request != section->end(); ++request) {
		if (!request.value().is_object()) {
			std::cout << "sample '" << request.key() << "' must be an object." << std::endl;
			return false;
		}
		static const std::vector<std::string> keys = {"type", "seed", "mol", "size"};
		for (auto item = request.value().begin(); item != request.value().end(); ++item) {
			if (ContainsValue(keys, item.key())) continue;
			std::cout << "sample property '" << item.key() << "' is unknown." << std::endl;
			return false;
		}
		try {
			if (request.value().value("type", std::string{}) != "xyz") {
				std::cout << "sample '" << request.key() << "' currently requires type 'xyz'." << std::endl;
				return false;
			}
			const std::string molecule_name = request.value().value("mol", std::string{});
			auto molecule = std::find_if(molecules.begin(), molecules.end(), [&](const auto& candidate) { return candidate->name == molecule_name; });
			if (molecule == molecules.end()) {
				std::cout << "sample '" << request.key() << "' references unknown molecule '" << molecule_name << "'." << std::endl;
				return false;
			}
			const int size = request.value().value("size", 1);
			if (size < 1) {
				std::cout << "sample '" << request.key() << "' requires size larger than zero." << std::endl;
				return false;
			}
			uint64_t seed;
			if (request.value().contains("seed")) seed = request.value().at("seed").get<uint64_t>();
			else {
				// The chosen value is written to the output, so an initially random request can be replayed exactly.
				std::random_device random;
				seed = (static_cast<uint64_t>(random()) << 32) ^ random();
			}
			samples.emplace_back(request.key(), molecule->get(), seed, size);
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid input for sample '" << request.key() << "': " << error.what() << std::endl;
			return false;
		}
	}
	return true;
}

}
