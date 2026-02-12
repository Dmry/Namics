#ifndef COMPLEX_PROPAGATORxH
#define COMPLEX_PROPAGATORxH

#include "../namics.h"
#include "../lattice.h"
#include "stl_typedef.h"

class ComplexPropagator {
public:
	ComplexPropagator(Lattice* lat);
	~ComplexPropagator();

	void compute_boltzmann(Real* W_R, Real* W_I, Real* G1_R, Real* G1_I);

	void propagate_step(Real* dst_R, Real* dst_I,
	                    Real* src_R, Real* src_I,
	                    Real* G1_R, Real* G1_I);

private:
	Lattice* lat;
	int M;
	stl::device_vector<Real> buf_R;
	stl::device_vector<Real> buf_I;
};

#endif
