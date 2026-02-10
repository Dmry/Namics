#ifndef CL_DYNAMICSxH
#define CL_DYNAMICSxH

#include "../namics.h"
#include "../lattice.h"
#include "../molecule.h"
#include "../segment.h"
#include "../system.h"
#include "complex_propagator.h"
#include "component.h"
#include <vector>
#include <memory>
#include <random>

using std::vector;
using std::unique_ptr;
using std::shared_ptr;

// Complex Langevin dynamics engine.
// Evolves complex chemical potential fields W = W_R + i*W_I via:
//   dW_R/dt = force_R + noise      (stochastic)
//   dW_I/dt = force_I              (deterministic)
// where force = chi*phi + xi - W (relaxation toward SCF saddle point).
class CL_Dynamics {
public:
	CL_Dynamics(Lattice* lat, const vector<Molecule*>& Mol, const vector<Segment*>& Seg,
	            System* Sys, Real dt, Real concentration, int seed);
	~CL_Dynamics();

	void initialize(vector<shared_ptr<IComponent>>& components);
	void step();

	Real* density_real(int comp) { return phi_R[comp]; }
	Real* density_imag(int comp) { return phi_I[comp]; }
	int n_comp() const { return n_comp_; }

private:
	Lattice* lat;
	vector<Molecule*> Mol;
	vector<Segment*> Seg;
	System* Sys;

	unique_ptr<ComplexPropagator> prop;

	int M;
	int n_comp_;
	int n_mol_;
	int n_seg_;
	Real dt_;
	Real noise_scale; // sqrt(2*dt/C)

	// Per-component (size M each)
	vector<Real*> W_R, W_I;
	vector<Real*> phi_R, phi_I;
	vector<Real*> force_R, force_I;

	// Per-segment-type Boltzmann factors (only allocated for component segments)
	vector<Real*> G1_R, G1_I;

	// Chain propagation work arrays
	Real* Gg_f_R;  // forward (max_chainlength * M)
	Real* Gg_f_I;
	Real* Gg_b_R;  // backward (2 * M, ping-pong)
	Real* Gg_b_I;
	int max_chainlength;

	// Incompressibility
	Real* pressure_R;
	Real* pressure_I;

	// Mappings
	vector<int> sysmolmon;    // comp idx -> seg idx
	vector<int> seg_to_comp;  // seg idx -> comp idx (-1 if frozen/absent)

	// Noise
	std::mt19937 rng;
	std::normal_distribution<Real> normal_dist;

	void compute_boltzmann_all();
	void compute_densities();
	void propagate_molecule(int mol_idx);
	void compute_forces();
	void update_fields();
};

#endif
