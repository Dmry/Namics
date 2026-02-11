#ifndef CL_DYNAMICSxH
#define CL_DYNAMICSxH

#include "../namics.h"
#include "../lattice.h"
#include "../molecule.h"
#include "../segment.h"
#include "../system.h"
#include "complex_propagator.h"
#include "component.h"
#include "stl_typedef.h"
#include <random>

class CL_Dynamics {
public:
	CL_Dynamics(Lattice* lat, const vector<Molecule*>& Mol, const vector<Segment*>& Seg,
	            System* Sys, Real dt, Real concentration, int seed);
	~CL_Dynamics();

	void initialize(vector<shared_ptr<IComponent>>& components);
	void step();

	Real* density_real(int comp);
	Real* density_imag(int comp);
	void sync_density_to_host();
	Real* density_real_host(int comp);
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
	Real C_;           // concentration parameter
	Real noise_scale;  // sqrt(2*dt/C)
	Real kappa_;       // compressibility penalty

	vector<stl::device_vector<Real>> W_R, W_I;       // per-component, size M each
	vector<stl::device_vector<Real>> phi_R, phi_I;
	vector<stl::device_vector<Real>> force_R, force_I;

	vector<stl::device_vector<Real>> G1_R, G1_I;     // Boltzmann factors per segment type

	// chain propagation
	stl::device_vector<Real> Gg_f_R, Gg_f_I;         // forward (max_chainlength * M)
	stl::device_vector<Real> Gg_b_R, Gg_b_I;         // backward (2 * M, ping-pong)
	int max_chainlength;

	stl::device_vector<Real> pressure_R, pressure_I;  // work buffer for total density

	vector<int> sysmolmon;        // comp idx -> seg idx
	vector<int> seg_to_comp;      // seg idx -> comp idx (-1 if absent)

	std::mt19937 rng;
	std::normal_distribution<Real> normal_dist;
	stl::host_vector<Real> h_noise;
	stl::device_vector<Real> d_noise;
	vector<stl::host_vector<Real>> h_phi_R;           // host copies for output

	void compute_boltzmann_all();
	void compute_densities();
	void propagate_molecule(int mol_idx);
	void compute_forces();
	void update_fields();
};

#endif
