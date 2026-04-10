#ifndef MOLECULExH
#define MOLECULExH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"

class Molecule {
public:
	struct Node {
		int segment = -1;
		int repeat = 1;
		std::vector<std::vector<Node>> branches;
	};
	struct SegmentOccurrence {
		int segment = -1;
		int segment_type_index = -1;
		std::vector<int> children;
	};
	struct OutputRequest {
		bool segment_density = false;
		bool ranked_density = false;
	};
	using Topology = std::vector<Node>;

	Molecule(Lattice*,std::span<const std::unique_ptr<Segment>>,std::string);

	std::string name;
	std::string composition;
	std::span<const std::unique_ptr<Segment>> Seg;
	Lattice* lat = nullptr;
	Topology topology; // parsed molecule graph
	std::vector<int> segment_types;
	Real Mu = 0;
	Real theta = 0;
	Real phibulk = 0;
	std::string freedom;
	Real n = 0;
	Real GN = 0;
	Real norm = 0;
	int chainlength = 0;
	std::vector<SegmentOccurrence> segment_path;
	std::vector<Real> mu_state;
	std::vector<Real> phi;
	std::vector<Real> phi_ranked;
	std::vector<Real> phitot;
	std::vector<Real> q_forward;
	std::vector<Real> G_unity;
	OutputRequest output_request;
	ParameterStore OUTPUT;
	std::span<const int> SegmentTypes() const noexcept { return segment_types; }
	template<typename F>
	void ForEachNode(F&& f) const { WalkNode(topology, f); }
	template<typename F>
	void ForEachOccurrence(F&& f) const { WalkOccurrence(topology, f); }
	void PushOutput();
	std::span<Real> GetPointer(int);

	bool IsPinned(void);
	bool IsCharged(void);
	Real Charge(void);
	void AllocateMemory(void);
	void PrepareForCalculations(std::span<const Real>);
	void ComputeGN();
	void FinalizeOutputs();
	Real fraction(int);
	void AccumulateDensity(std::span<Real> system_phitot);

private:
	template<typename F>
	static void WalkNode(const Topology& chain, F&& f) {
		for (const auto& node : chain) {
			f(node);
			for (const auto& branch : node.branches) WalkNode(branch, f);
		}
	}
	template<typename F>
	static void WalkOccurrence(const Topology& chain, F&& f) {
		for (const auto& node : chain) {
			for (int repeat = 0; repeat < node.repeat; ++repeat) f(node);
			for (const auto& branch : node.branches) WalkOccurrence(branch, f);
		}
	}

	void AddSegmentDensity(int, std::span<const Real>, bool, std::span<Real>, std::span<Real>);
	void PropagateForward();
	void PropagateBackward(int, std::span<const Real>, bool, std::span<Real>, std::span<Real>);
	void AccumulateDensity(bool, std::span<Real>, std::span<Real>);
};

namespace molecule_factory {
std::unique_ptr<Molecule> CreateChecked(const Input&, Lattice*, std::span<const std::unique_ptr<Segment>>, const std::string&, int);
}

#endif
