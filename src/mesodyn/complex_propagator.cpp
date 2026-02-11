#include "complex_propagator.h"
#include <cmath>

struct boltzmann_functor {
	Real *W_R, *W_I, *G1_R, *G1_I;
	boltzmann_functor(Real* wr, Real* wi, Real* gr, Real* gi)
		: W_R(wr), W_I(wi), G1_R(gr), G1_I(gi) {}
	DEVICE_LAMBDA void operator()(int i) const {
		Real e = exp(-W_R[i]);
		G1_R[i] =  e * cos(W_I[i]);
		G1_I[i] = -e * sin(W_I[i]);
	}
};

struct complex_multiply_functor {
	Real *dst_R, *dst_I, *a_R, *a_I, *b_R, *b_I;
	complex_multiply_functor(Real* dr, Real* di, Real* ar, Real* ai, Real* br, Real* bi)
		: dst_R(dr), dst_I(di), a_R(ar), a_I(ai), b_R(br), b_I(bi) {}
	DEVICE_LAMBDA void operator()(int i) const {
		dst_R[i] = a_R[i] * b_R[i] - a_I[i] * b_I[i];
		dst_I[i] = a_R[i] * b_I[i] + a_I[i] * b_R[i];
	}
};

ComplexPropagator::ComplexPropagator(Lattice* lat_)
	: lat(lat_), M(lat_->M), buf_R(M, 0.0), buf_I(M, 0.0) {}

ComplexPropagator::~ComplexPropagator() {}

void ComplexPropagator::compute_boltzmann(Real* W_R, Real* W_I,
                                          Real* G1_R, Real* G1_I) {
	parallel_for(M, boltzmann_functor(W_R, W_I, G1_R, G1_I));
}

void ComplexPropagator::propagate_step(Real* dst_R, Real* dst_I,
                                       Real* src_R, Real* src_I,
                                       Real* G1_R, Real* G1_I) {
	lat->set_bounds(src_R);
	lat->set_bounds(src_I);

	Real* br = raw_ptr(buf_R);
	Real* bi = raw_ptr(buf_I);
	lat->Side(br, src_R, M);
	lat->Side(bi, src_I, M);

	parallel_for(M, complex_multiply_functor(dst_R, dst_I, G1_R, G1_I, br, bi));
}
