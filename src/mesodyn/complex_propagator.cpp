#include "complex_propagator.h"
#include "../tools.h"
#include <cmath>

ComplexPropagator::ComplexPropagator(Lattice* lat_) : lat(lat_), M(lat_->M) {
	buf_R = (Real*)malloc(M * sizeof(Real));
	buf_I = (Real*)malloc(M * sizeof(Real));
}

ComplexPropagator::~ComplexPropagator() {
	free(buf_R);
	free(buf_I);
}

void ComplexPropagator::compute_boltzmann(Real* W_R, Real* W_I,
                                          Real* G1_R, Real* G1_I) {
	for (int i = 0; i < M; i++) {
		Real e = exp(-W_R[i]);
		G1_R[i] =  e * cos(W_I[i]);
		G1_I[i] = -e * sin(W_I[i]);
	}
}

void ComplexPropagator::propagate_step(Real* dst_R, Real* dst_I,
                                       Real* src_R, Real* src_I,
                                       Real* G1_R, Real* G1_I) {
	lat->set_bounds(src_R);
	lat->set_bounds(src_I);

	lat->Side(buf_R, src_R, M);
	lat->Side(buf_I, src_I, M);

	for (int i = 0; i < M; i++) {
		dst_R[i] = G1_R[i] * buf_R[i] - G1_I[i] * buf_I[i];
		dst_I[i] = G1_R[i] * buf_I[i] + G1_I[i] * buf_R[i];
	}
}
