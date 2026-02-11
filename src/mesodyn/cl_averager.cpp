#include "cl_averager.h"
#include <cmath>
#include <algorithm>

CL_Averager::CL_Averager(int n_comp, int M, size_t equilibration)
	: n_comp_(n_comp), M_(M), equilibration_(equilibration), count_(0)
{
	for (int i = 0; i < n_comp; i++) {
		sum_.emplace_back(M, 0.0);
		sum_sq_.emplace_back(M, 0.0);
		mean_.emplace_back(M, 0.0);
		var_.emplace_back(M, 0.0);
	}
}

CL_Averager::~CL_Averager() {}

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
		std::fill(sum_[i].begin(), sum_[i].end(), 0.0);
		std::fill(sum_sq_[i].begin(), sum_sq_[i].end(), 0.0);
		std::fill(mean_[i].begin(), mean_[i].end(), 0.0);
		std::fill(var_[i].begin(), var_[i].end(), 0.0);
	}
}
