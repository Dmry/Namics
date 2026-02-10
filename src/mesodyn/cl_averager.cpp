#include "cl_averager.h"
#include "../tools.h"
#include <cmath>

CL_Averager::CL_Averager(int n_comp, int M, size_t equilibration)
	: n_comp_(n_comp), M_(M), equilibration_(equilibration), count_(0)
{
	sum_.resize(n_comp);
	sum_sq_.resize(n_comp);
	mean_.resize(n_comp);
	var_.resize(n_comp);

	for (int i = 0; i < n_comp; i++) {
		sum_[i]    = (Real*)calloc(M, sizeof(Real));
		sum_sq_[i] = (Real*)calloc(M, sizeof(Real));
		mean_[i]   = (Real*)calloc(M, sizeof(Real));
		var_[i]    = (Real*)calloc(M, sizeof(Real));
	}
}

CL_Averager::~CL_Averager() {
	for (int i = 0; i < n_comp_; i++) {
		free(sum_[i]);
		free(sum_sq_[i]);
		free(mean_[i]);
		free(var_[i]);
	}
}

void CL_Averager::accumulate(vector<Real*>& phi_R, size_t t) {
	if (t < equilibration_) return;

	count_++;
	Real inv_n = 1.0 / count_;

	for (int i = 0; i < n_comp_; i++) {
		for (int r = 0; r < M_; r++) {
			Real val = phi_R[i][r];
			sum_[i][r]    += val;
			sum_sq_[i][r] += val * val;
			mean_[i][r] = sum_[i][r] * inv_n;
			var_[i][r]  = sum_sq_[i][r] * inv_n - mean_[i][r] * mean_[i][r];
		}
	}
}

void CL_Averager::reset() {
	count_ = 0;
	for (int i = 0; i < n_comp_; i++) {
		Zero(sum_[i], M_);
		Zero(sum_sq_[i], M_);
		Zero(mean_[i], M_);
		Zero(var_[i], M_);
	}
}
