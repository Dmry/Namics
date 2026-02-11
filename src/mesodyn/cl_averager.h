#ifndef CL_AVERAGERxH
#define CL_AVERAGERxH

#include "../namics.h"
#include <vector>

class CL_Averager {
public:
	CL_Averager(int n_comp, int M, size_t equilibration);
	~CL_Averager();

	void accumulate(vector<Real*>& phi_R, size_t t);

	Real* mean(int comp) { return mean_[comp].data(); }
	Real* variance(int comp) { return var_[comp].data(); }

	size_t sample_count() const { return count_; }

	void reset();

private:
	int n_comp_;
	int M_;
	size_t equilibration_;
	size_t count_;

	vector<vector<Real>> sum_;
	vector<vector<Real>> sum_sq_;
	vector<vector<Real>> mean_;
	vector<vector<Real>> var_;
};

#endif
