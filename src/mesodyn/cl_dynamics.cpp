#include "cl_dynamics.h"
#include "../tools.h"
#include <cmath>
#include <algorithm>
#include <iostream>

CL_Dynamics::CL_Dynamics(Lattice* lat_,
                         const vector<Molecule*>& Mol_,
                         const vector<Segment*>& Seg_,
                         System* Sys_,
                         Real dt,
                         Real concentration,
                         int seed)
	: lat(lat_), Mol(Mol_), Seg(Seg_), Sys(Sys_),
	  M(lat_->M), dt_(dt),
	  noise_scale(sqrt(2.0 * dt / concentration)),
	  rng(seed), normal_dist(0.0, 1.0)
{
	prop = unique_ptr<ComplexPropagator>(new ComplexPropagator(lat));

	n_comp_ = Sys->SysMolMonList.size();
	n_mol_  = Mol.size();
	n_seg_  = Seg.size();

	// Build segment <-> component mappings
	sysmolmon.resize(n_comp_);
	for (int i = 0; i < n_comp_; i++)
		sysmolmon[i] = Sys->SysMolMonList[i];

	seg_to_comp.assign(n_seg_, -1);
	for (int i = 0; i < n_comp_; i++)
		seg_to_comp[sysmolmon[i]] = i;

	// Allocate per-component arrays
	W_R.resize(n_comp_);     W_I.resize(n_comp_);
	phi_R.resize(n_comp_);   phi_I.resize(n_comp_);
	force_R.resize(n_comp_); force_I.resize(n_comp_);

	for (int i = 0; i < n_comp_; i++) {
		W_R[i]     = (Real*)calloc(M, sizeof(Real));
		W_I[i]     = (Real*)calloc(M, sizeof(Real));
		phi_R[i]   = (Real*)calloc(M, sizeof(Real));
		phi_I[i]   = (Real*)calloc(M, sizeof(Real));
		force_R[i] = (Real*)calloc(M, sizeof(Real));
		force_I[i] = (Real*)calloc(M, sizeof(Real));
	}

	// Boltzmann factors indexed by segment type (only allocate for components)
	G1_R.assign(n_seg_, nullptr);
	G1_I.assign(n_seg_, nullptr);
	for (int i = 0; i < n_comp_; i++) {
		int seg = sysmolmon[i];
		G1_R[seg] = (Real*)calloc(M, sizeof(Real));
		G1_I[seg] = (Real*)calloc(M, sizeof(Real));
	}

	pressure_R = (Real*)calloc(M, sizeof(Real));
	pressure_I = (Real*)calloc(M, sizeof(Real));

	// Forward propagator sized for longest chain
	max_chainlength = 1;
	for (int j = 0; j < n_mol_; j++)
		if (Mol[j]->chainlength > max_chainlength)
			max_chainlength = Mol[j]->chainlength;

	Gg_f_R = (Real*)calloc(max_chainlength * M, sizeof(Real));
	Gg_f_I = (Real*)calloc(max_chainlength * M, sizeof(Real));
	Gg_b_R = (Real*)calloc(2 * M, sizeof(Real));
	Gg_b_I = (Real*)calloc(2 * M, sizeof(Real));
}

CL_Dynamics::~CL_Dynamics() {
	for (int i = 0; i < n_comp_; i++) {
		free(W_R[i]);     free(W_I[i]);
		free(phi_R[i]);   free(phi_I[i]);
		free(force_R[i]); free(force_I[i]);
	}
	for (int i = 0; i < n_comp_; i++) {
		int seg = sysmolmon[i];
		free(G1_R[seg]); free(G1_I[seg]);
	}
	free(pressure_R); free(pressure_I);
	free(Gg_f_R); free(Gg_f_I);
	free(Gg_b_R); free(Gg_b_I);
}

void CL_Dynamics::initialize(vector<shared_ptr<IComponent>>& components) {
	// Initialize W_R from mean-field relation: W_i = Sigma_k chi_{ik} * phi_k
	for (int i = 0; i < n_comp_; i++) {
		Zero(W_R[i], M);
		Zero(W_I[i], M);

		int seg_i = sysmolmon[i];
		for (int k = 0; k < n_seg_; k++) {
			Real chi = Sys->CHI[seg_i * n_seg_ + k];
			if (chi != 0.0) {
				int comp_k = seg_to_comp[k];
				if (comp_k >= 0) {
					Real* rho = (Real*)components[comp_k]->rho;
					YplusisCtimesX(W_R[i], rho, chi, M);
				} else {
					YplusisCtimesX(W_R[i], Seg[k]->phi_side, chi, M);
				}
			}
		}
	}
}

void CL_Dynamics::compute_boltzmann_all() {
	for (int i = 0; i < n_comp_; i++) {
		int seg = sysmolmon[i];
		prop->compute_boltzmann(W_R[i], W_I[i], G1_R[seg], G1_I[seg]);
	}
}

void CL_Dynamics::propagate_molecule(int mol_idx) {
	Molecule* mol = Mol[mol_idx];
	int chainlength = mol->chainlength;
	int n_blocks = mol->n_mon.size();

	// Forward propagation
	int s = 0;
	for (int b = 0; b < n_blocks; b++) {
		int seg = mol->mon_nr[b];
		int N = mol->n_mon[b];

		for (int k = 0; k < N; k++) {
			if (s == 0) {
				// Initialize: Gg_f(0) = G1
				Cp(Gg_f_R, G1_R[seg], M);
				Cp(Gg_f_I, G1_I[seg], M);
			} else {
				prop->propagate_step(
					Gg_f_R + s * M, Gg_f_I + s * M,
					Gg_f_R + (s - 1) * M, Gg_f_I + (s - 1) * M,
					G1_R[seg], G1_I[seg]);
			}
			s++;
		}
	}

	// Partition function: Q = WeightedSum(Gg_f at last segment)
	Real Q_R = lat->WeightedSum(Gg_f_R + (chainlength - 1) * M);
	Real Q_I = lat->WeightedSum(Gg_f_I + (chainlength - 1) * M);

	// norm = phibulk / (chainlength * Q)  [complex division: real / complex]
	Real phibulk = mol->phibulk;
	Real a = phibulk / chainlength;
	Real Q_mag2 = Q_R * Q_R + Q_I * Q_I;
	if (Q_mag2 < 1e-30) Q_mag2 = 1e-30;
	Real norm_R =  a * Q_R / Q_mag2;
	Real norm_I = -a * Q_I / Q_mag2;

	// Backward propagation + density accumulation
	s = chainlength - 1;
	for (int b = n_blocks - 1; b >= 0; b--) {
		int seg = mol->mon_nr[b];
		int N = mol->n_mon[b];
		int comp = seg_to_comp[seg];

		for (int k = N - 1; k >= 0; k--) {
			if (s == chainlength - 1) {
				Cp(Gg_b_R + (s % 2) * M, G1_R[seg], M);
				Cp(Gg_b_I + (s % 2) * M, G1_I[seg], M);
			} else {
				prop->propagate_step(
					Gg_b_R + (s % 2) * M, Gg_b_I + (s % 2) * M,
					Gg_b_R + ((s + 1) % 2) * M, Gg_b_I + ((s + 1) % 2) * M,
					G1_R[seg], G1_I[seg]);
			}

			// phi += norm * Gg_f(s) * Gg_b(s) / G1
			if (comp >= 0) {
				Real* f_R = Gg_f_R + s * M;
				Real* f_I = Gg_f_I + s * M;
				Real* b_R = Gg_b_R + (s % 2) * M;
				Real* b_I = Gg_b_I + (s % 2) * M;

				for (int r = 0; r < M; r++) {
					// product = Gg_f * Gg_b
					Real p_R = f_R[r] * b_R[r] - f_I[r] * b_I[r];
					Real p_I = f_R[r] * b_I[r] + f_I[r] * b_R[r];

					// product / G1
					Real g_R = G1_R[seg][r];
					Real g_I = G1_I[seg][r];
					Real g2 = g_R * g_R + g_I * g_I;
					if (g2 < 1e-30) g2 = 1e-30;
					Real d_R = (p_R * g_R + p_I * g_I) / g2;
					Real d_I = (p_I * g_R - p_R * g_I) / g2;

					// phi += norm * d
					phi_R[comp][r] += norm_R * d_R - norm_I * d_I;
					phi_I[comp][r] += norm_R * d_I + norm_I * d_R;
				}
			}

			s--;
		}
	}
}

void CL_Dynamics::compute_densities() {
	for (int i = 0; i < n_comp_; i++) {
		Zero(phi_R[i], M);
		Zero(phi_I[i], M);
	}
	for (int j = 0; j < n_mol_; j++)
		propagate_molecule(j);
}

void CL_Dynamics::compute_forces() {
	// force_i = Sigma_k chi_{ik} * phi_k + xi - W_i
	// xi = (Sigma_i W_i - Sigma_i (chi*phi)_i) / n_comp

	// First: accumulate chi*phi for each component
	for (int i = 0; i < n_comp_; i++) {
		Zero(force_R[i], M);
		Zero(force_I[i], M);

		int seg_i = sysmolmon[i];
		for (int k = 0; k < n_seg_; k++) {
			Real chi = Sys->CHI[seg_i * n_seg_ + k];
			if (chi != 0.0) {
				int comp_k = seg_to_comp[k];
				if (comp_k >= 0) {
					YplusisCtimesX(force_R[i], phi_R[comp_k], chi, M);
					YplusisCtimesX(force_I[i], phi_I[comp_k], chi, M);
				} else {
					// Frozen segment: real density only
					YplusisCtimesX(force_R[i], Seg[k]->phi_side, chi, M);
				}
			}
		}
	}

	// Pressure: xi = (Sigma_i W_i - Sigma_i chi*phi_i) / n_comp
	Zero(pressure_R, M);
	Zero(pressure_I, M);
	for (int i = 0; i < n_comp_; i++) {
		for (int r = 0; r < M; r++) {
			pressure_R[r] += W_R[i][r] - force_R[i][r];
			pressure_I[r] += W_I[i][r] - force_I[i][r];
		}
	}
	Norm(pressure_R, 1.0 / n_comp_, M);
	Norm(pressure_I, 1.0 / n_comp_, M);

	// force_i = chi*phi + xi - W_i
	for (int i = 0; i < n_comp_; i++) {
		for (int r = 0; r < M; r++) {
			force_R[i][r] += pressure_R[r] - W_R[i][r];
			force_I[i][r] += pressure_I[r] - W_I[i][r];
		}
	}
}

void CL_Dynamics::update_fields() {
	for (int i = 0; i < n_comp_; i++) {
		for (int r = 0; r < M; r++) {
			Real eta = normal_dist(rng);
			W_R[i][r] += dt_ * force_R[i][r] + noise_scale * eta;
			W_I[i][r] += dt_ * force_I[i][r];
		}
	}
}

void CL_Dynamics::step() {
	compute_boltzmann_all();
	compute_densities();
	compute_forces();
	update_fields();
}
