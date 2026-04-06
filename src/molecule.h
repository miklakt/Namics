#ifndef MOLECULExH
#define MOLECULExH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"
#include "tools.h"

class Molecule {
public:
	struct Node {
		int segment = -1;
		int repeat = 1;
		std::vector<std::vector<Node>> branches;
	};
	using Topology = std::vector<Node>;

	Molecule(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::string);
	virtual ~Molecule();

	std::string name;
	bool all_molecule;
	const Input* In;
	std::span<const std::unique_ptr<Segment>> Seg;
	Lattice* lat;
	Topology topology; // parsed molecule graph
	std::vector<int> MolMonList;
	int start;
	Real Mu;
	Real theta;
	Real phibulk;
	std::string freedom;
	Real n;
	Real GN;
	Real norm;
	int chainlength,N;
	std::vector<int> Gnr; //generation-number
	std::vector<int> first_s;
	std::vector<int> last_s;
	std::vector<int> first_b;
	std::vector<int> last_b;
	std::vector<int> mon_nr;
	std::vector<int> n_mon;
	std::vector<int> molmon_nr;
	std::vector<Real> mu_state;
	std::vector<Real> phi;
	std::vector<Real> phi_ranked;
	std::vector<Real> phitot;
	std::vector<Real> Gg_f;
	std::vector<Real> Gg_b;
	std::vector<Real> UNITY;
	Real B;
	ParameterStore OUTPUT;
	std::span<const int> SegmentIndices() const noexcept { return MolMonList; }
	template<typename F>
	void ForEachNode(F&& f) const { WalkNode(topology, 0, f); }
	void PushOutput();
	std::span<Real> GetPointer(int);

	bool IsPinned(void);
	bool IsCharged(void);
	Real Charge(void);
	bool HasOutputProperty(const std::string&) const;
	void DeAllocateMemory(void);
	void AllocateMemory(void);
	bool PrepareForCalculations(std::span<const Real>);
	virtual bool ComputePhi();
	virtual void FinalizeOutputs();
	virtual Real fraction(int);

	Real* propagate_forward(Real*,int&,int,int,int);
	void propagate_backward(Real*,int&,int,int,std::span<Real> ranked_phi = {});

private:
	template<typename F>
	static void WalkNode(const Topology& chain, int depth, F&& f) {
		for (const auto& node : chain) {
			f(node, depth);
			for (const auto& branch : node.branches) WalkNode(branch, depth + 1, f);
		}
	}

	Real* ForwardBranch(int generation, int &s);
	void BackwardBranch(int generation, int &s, std::span<Real> ranked_phi = {});
	bool ComputePhiRanked(std::span<Real> ranked_phi);
};

namespace molecule_factory {
std::unique_ptr<Molecule> CreateChecked(const Input&, Lattice*, std::span<const std::unique_ptr<Segment>>, const std::string&, int);
}

#endif
