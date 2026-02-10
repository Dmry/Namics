#ifndef COMPLEX_PROPAGATORxH
#define COMPLEX_PROPAGATORxH

#include "../namics.h"
#include "../lattice.h"

// Complex Scheutjens-Fleer lattice propagation using paired Real arrays.
// Reads lattice geometry (lambda weights, strides) from Lat via Side().
class ComplexPropagator {
public:
	ComplexPropagator(Lattice* lat);
	~ComplexPropagator();

	// G1 = exp(-(W_R + i*W_I))
	void compute_boltzmann(Real* W_R, Real* W_I, Real* G1_R, Real* G1_I);

	// One propagation step: dst = G1 * Side(src)
	// src arrays are modified by set_bounds (boundary fill).
	void propagate_step(Real* dst_R, Real* dst_I,
	                    Real* src_R, Real* src_I,
	                    Real* G1_R, Real* G1_I);

private:
	Lattice* lat;
	int M;
	Real* buf_R;
	Real* buf_I;
};

#endif
