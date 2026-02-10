#ifndef CL_AVERAGERxH
#define CL_AVERAGERxH

#include "../namics.h"
#include <vector>

using std::vector;

// Running time-average of CL density fields (real part).
// Accumulates after a configurable equilibration period.
// Tracks mean and variance for error estimation.
class CL_Averager {
public:
	CL_Averager(int n_comp, int M, size_t equilibration);
	~CL_Averager();

	// Accumulate one sample (real part of density for each component)
	void accumulate(vector<Real*>& phi_R, size_t t);

	// Access averaged density for component i
	Real* mean(int comp) { return mean_[comp]; }

	// Access variance for component i
	Real* variance(int comp) { return var_[comp]; }

	size_t sample_count() const { return count_; }

	void reset();

private:
	int n_comp_;
	int M_;
	size_t equilibration_;
	size_t count_;

	vector<Real*> sum_;     // running sum
	vector<Real*> sum_sq_;  // running sum of squares
	vector<Real*> mean_;    // current mean
	vector<Real*> var_;     // current variance
};

#endif
