#include "cl_dynamics.h"
#include "../tools.h"
#include <cmath>
#include <algorithm>

#ifdef PAR_MESODYN_THRUST
#include <thrust/transform.h>
#endif

struct phi_accumulate_functor {
	Real *phi_R, *phi_I;
	Real *f_R, *f_I, *b_R, *b_I;
	Real *g_R, *g_I;
	Real norm_R, norm_I;
	phi_accumulate_functor(Real* pr, Real* pi, Real* fr, Real* fi,
	    Real* br, Real* bi, Real* gr, Real* gi, Real nr, Real ni)
		: phi_R(pr), phi_I(pi), f_R(fr), f_I(fi),
		  b_R(br), b_I(bi), g_R(gr), g_I(gi), norm_R(nr), norm_I(ni) {}
	DEVICE_LAMBDA void operator()(int r) const {
		Real p_R = f_R[r] * b_R[r] - f_I[r] * b_I[r];
		Real p_I = f_R[r] * b_I[r] + f_I[r] * b_R[r];
		Real g2 = g_R[r] * g_R[r] + g_I[r] * g_I[r];
		if (g2 < 1e-30) g2 = 1e-30;
		Real d_R = (p_R * g_R[r] + p_I * g_I[r]) / g2;
		Real d_I = (p_I * g_R[r] - p_R * g_I[r]) / g2;
		phi_R[r] += norm_R * d_R - norm_I * d_I;
		phi_I[r] += norm_R * d_I + norm_I * d_R;
	}
};

CL_Dynamics::CL_Dynamics(Lattice* lat_,
                         const vector<Molecule*>& Mol_,
                         const vector<Segment*>& Seg_,
                         System* Sys_,
                         Real dt,
                         Real concentration,
                         int seed)
	: lat(lat_), Mol(Mol_), Seg(Seg_), Sys(Sys_),
	  M(lat_->M), dt_(dt), C_(concentration),
	  noise_scale(sqrt(2.0 * dt / concentration)),
	  kappa_(100.0),
	  rng(seed), normal_dist(0.0, 1.0)
{
	prop = unique_ptr<ComplexPropagator>(new ComplexPropagator(lat));

	n_comp_ = Sys->SysMolMonList.size();
	n_mol_  = Mol.size();
	n_seg_  = Seg.size();

	sysmolmon.resize(n_comp_);
	for (int i = 0; i < n_comp_; i++)
		sysmolmon[i] = Sys->SysMolMonList[i];

	seg_to_comp.assign(n_seg_, -1);
	for (int i = 0; i < n_comp_; i++)
		seg_to_comp[sysmolmon[i]] = i;

	for (int i = 0; i < n_comp_; i++) {
		W_R.emplace_back(M, 0.0);     W_I.emplace_back(M, 0.0);
		phi_R.emplace_back(M, 0.0);   phi_I.emplace_back(M, 0.0);
		force_R.emplace_back(M, 0.0); force_I.emplace_back(M, 0.0);
	}

	G1_R.resize(n_seg_);
	G1_I.resize(n_seg_);
	for (int i = 0; i < n_comp_; i++) {
		int seg = sysmolmon[i];
		G1_R[seg].resize(M, 0.0);
		G1_I[seg].resize(M, 0.0);
	}

	pressure_R.resize(M, 0.0);
	pressure_I.resize(M, 0.0);

	max_chainlength = 1;
	for (int j = 0; j < n_mol_; j++)
		if (Mol[j]->chainlength > max_chainlength)
			max_chainlength = Mol[j]->chainlength;

	Gg_f_R.resize(max_chainlength * M, 0.0);
	Gg_f_I.resize(max_chainlength * M, 0.0);
	Gg_b_R.resize(2 * M, 0.0);
	Gg_b_I.resize(2 * M, 0.0);

	d_noise.resize(M);
	h_noise.resize(M);
	for (int i = 0; i < n_comp_; i++)
		h_phi_R.emplace_back(M, 0.0);
}

CL_Dynamics::~CL_Dynamics() {}

void CL_Dynamics::initialize(vector<shared_ptr<IComponent>>& components) {
	Real total_phibulk = 0;
	int solvent_mol = -1;
	for (int j = 0; j < n_mol_; j++) {
		if (Mol[j]->phibulk > 0) {
			total_phibulk += Mol[j]->phibulk;
		} else if (Mol[j]->freedom == "solvent") {
			solvent_mol = j;
		} else if (Mol[j]->theta > 0) {
			Mol[j]->phibulk = Mol[j]->theta / lat->volume;
			total_phibulk += Mol[j]->phibulk;
		}
	}
	if (solvent_mol >= 0)
		Mol[solvent_mol]->phibulk = 1.0 - total_phibulk;

	for (int i = 0; i < n_comp_; i++) {
		stl::fill(EXEC_PAR W_R[i].begin(), W_R[i].end(), 0.0);
		stl::fill(EXEC_PAR W_I[i].begin(), W_I[i].end(), 0.0);

		int seg_i = sysmolmon[i];
		for (int k = 0; k < n_seg_; k++) {
			Real chi = Sys->CHI[seg_i * n_seg_ + k];
			if (chi != 0.0) {
				int comp_k = seg_to_comp[k];
				if (comp_k >= 0) {
					stl::transform(EXEC_PAR
						components[comp_k]->rho.begin(),
						components[comp_k]->rho.end(),
						W_R[i].begin(), W_R[i].begin(),
						saxpy_functor(chi));
				} else {
					YplusisCtimesX(raw_ptr(W_R[i]), Seg[k]->phi_side, chi, M);
				}
			}
		}
	}
}

void CL_Dynamics::compute_boltzmann_all() {
	for (int i = 0; i < n_comp_; i++) {
		int seg = sysmolmon[i];
		prop->compute_boltzmann(raw_ptr(W_R[i]), raw_ptr(W_I[i]),
		                        raw_ptr(G1_R[seg]), raw_ptr(G1_I[seg]));
	}
}

void CL_Dynamics::propagate_molecule(int mol_idx) {
	Molecule* mol = Mol[mol_idx];
	int chainlength = mol->chainlength;
	int n_blocks = mol->n_mon.size();

	Real* gf_R = raw_ptr(Gg_f_R);
	Real* gf_I = raw_ptr(Gg_f_I);
	Real* gb_R = raw_ptr(Gg_b_R);
	Real* gb_I = raw_ptr(Gg_b_I);

	// forward propagation
	int s = 0;
	for (int b = 0; b < n_blocks; b++) {
		int seg = mol->mon_nr[b];
		int N = mol->n_mon[b];
		Real* g1r = raw_ptr(G1_R[seg]);
		Real* g1i = raw_ptr(G1_I[seg]);

		for (int k = 0; k < N; k++) {
			if (s == 0) {
				stl::copy(EXEC_PAR G1_R[seg].begin(), G1_R[seg].end(),
					Gg_f_R.begin());
				stl::copy(EXEC_PAR G1_I[seg].begin(), G1_I[seg].end(),
					Gg_f_I.begin());
			} else {
				prop->propagate_step(
					gf_R + s * M, gf_I + s * M,
					gf_R + (s - 1) * M, gf_I + (s - 1) * M,
					g1r, g1i);
			}
			s++;
		}
	}

	Real Q_R = lat->WeightedSum(gf_R + (chainlength - 1) * M);
	lat->set_bounds(gf_R + (chainlength - 1) * M);
	Real Q_I = lat->WeightedSum(gf_I + (chainlength - 1) * M);
	lat->set_bounds(gf_I + (chainlength - 1) * M);

	Real phibulk = mol->phibulk;
	Real a = phibulk * lat->volume / chainlength;
	Real Q_mag2 = Q_R * Q_R + Q_I * Q_I;
	if (Q_mag2 < 1e-30) Q_mag2 = 1e-30;
	Real norm_R =  a * Q_R / Q_mag2;
	Real norm_I = -a * Q_I / Q_mag2;

	// backward pass + density accumulation
	s = chainlength - 1;
	for (int b = n_blocks - 1; b >= 0; b--) {
		int seg = mol->mon_nr[b];
		int N = mol->n_mon[b];
		int comp = seg_to_comp[seg];
		Real* g1r = raw_ptr(G1_R[seg]);
		Real* g1i = raw_ptr(G1_I[seg]);

		for (int k = N - 1; k >= 0; k--) {
			if (s == chainlength - 1) {
				stl::copy(EXEC_PAR G1_R[seg].begin(), G1_R[seg].end(),
					Gg_b_R.begin() + (s % 2) * M);
				stl::copy(EXEC_PAR G1_I[seg].begin(), G1_I[seg].end(),
					Gg_b_I.begin() + (s % 2) * M);
			} else {
				prop->propagate_step(
					gb_R + (s % 2) * M, gb_I + (s % 2) * M,
					gb_R + ((s + 1) % 2) * M, gb_I + ((s + 1) % 2) * M,
					g1r, g1i);
			}

			if (comp >= 0) {
				parallel_for(M, phi_accumulate_functor(
					raw_ptr(phi_R[comp]), raw_ptr(phi_I[comp]),
					gf_R + s * M, gf_I + s * M,
					gb_R + (s % 2) * M, gb_I + (s % 2) * M,
					g1r, g1i, norm_R, norm_I));
			}

			s--;
		}
	}
}

void CL_Dynamics::compute_densities() {
	for (int i = 0; i < n_comp_; i++) {
		stl::fill(EXEC_PAR phi_R[i].begin(), phi_R[i].end(), 0.0);
		stl::fill(EXEC_PAR phi_I[i].begin(), phi_I[i].end(), 0.0);
	}
	for (int j = 0; j < n_mol_; j++)
		propagate_molecule(j);
}

void CL_Dynamics::compute_forces() {
	for (int i = 0; i < n_comp_; i++) {
		stl::fill(EXEC_PAR force_R[i].begin(), force_R[i].end(), 0.0);
		stl::fill(EXEC_PAR force_I[i].begin(), force_I[i].end(), 0.0);

		int seg_i = sysmolmon[i];
		for (int k = 0; k < n_seg_; k++) {
			Real chi = Sys->CHI[seg_i * n_seg_ + k];
			if (chi != 0.0) {
				int comp_k = seg_to_comp[k];
				if (comp_k >= 0) {
					stl::transform(EXEC_PAR
						phi_R[comp_k].begin(), phi_R[comp_k].end(),
						force_R[i].begin(), force_R[i].begin(),
						saxpy_functor(chi));
					stl::transform(EXEC_PAR
						phi_I[comp_k].begin(), phi_I[comp_k].end(),
						force_I[i].begin(), force_I[i].begin(),
						saxpy_functor(chi));
				} else {
					YplusisCtimesX(raw_ptr(force_R[i]), Seg[k]->phi_side, chi, M);
				}
			}
		}
	}

	// total density
	stl::fill(EXEC_PAR pressure_R.begin(), pressure_R.end(), 0.0);
	stl::fill(EXEC_PAR pressure_I.begin(), pressure_I.end(), 0.0);
	for (int i = 0; i < n_comp_; i++) {
		stl::transform(EXEC_PAR phi_R[i].begin(), phi_R[i].end(), pressure_R.begin(), pressure_R.begin(), stl::plus<Real>());
		stl::transform(EXEC_PAR phi_I[i].begin(), phi_I[i].end(), pressure_I.begin(), pressure_I.begin(), stl::plus<Real>());
	}

	// compressibility penalty + field relaxation
	for (int i = 0; i < n_comp_; i++) {
		stl::transform(EXEC_PAR pressure_R.begin(), pressure_R.end(), force_R[i].begin(), force_R[i].begin(), compressibility_functor(kappa_));
		stl::transform(EXEC_PAR pressure_I.begin(), pressure_I.end(), force_I[i].begin(), force_I[i].begin(), saxpy_functor(kappa_));
		stl::transform(EXEC_PAR force_R[i].begin(), force_R[i].end(), W_R[i].begin(), force_R[i].begin(), stl::minus<Real>());
		stl::transform(EXEC_PAR force_I[i].begin(), force_I[i].end(), W_I[i].begin(), force_I[i].begin(), stl::minus<Real>());
	}
}

void CL_Dynamics::update_fields() {
	for (int i = 0; i < n_comp_; i++) {
		stl::transform(EXEC_PAR force_R[i].begin(), force_R[i].end(), W_R[i].begin(), W_R[i].begin(), saxpy_functor(dt_));
		stl::transform(EXEC_PAR force_I[i].begin(), force_I[i].end(), W_I[i].begin(), W_I[i].begin(), saxpy_functor(dt_));

		for (int r = 0; r < M; r++)
			h_noise[r] = noise_scale * normal_dist(rng);
		d_noise = h_noise;
		stl::transform(EXEC_PAR d_noise.begin(), d_noise.end(), W_R[i].begin(), W_R[i].begin(), stl::plus<Real>());
	}
}

Real* CL_Dynamics::density_real(int comp) { return raw_ptr(phi_R[comp]); }
Real* CL_Dynamics::density_imag(int comp) { return raw_ptr(phi_I[comp]); }

void CL_Dynamics::sync_density_to_host() {
	for (int i = 0; i < n_comp_; i++)
		h_phi_R[i] = phi_R[i];
}

Real* CL_Dynamics::density_real_host(int comp) {
	return h_phi_R[comp].data();
}

void CL_Dynamics::step() {
	compute_boltzmann_all();
	compute_densities();
	compute_forces();
	update_fields();
}
