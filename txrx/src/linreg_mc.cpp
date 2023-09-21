
#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <algorithm>
using namespace Eigen;
using namespace std;

//accelerate the part of cho_solve cho_factor of linreg_mc "linreg.py"
//sizes: rx_view = [N,  k], where bn = N - m + 1, bn + m - 1 = N;
//returns: scale_all [i- 1, 

extern"C" void ldlt_linreg(const float* _rx_view, const float* _sxy_seg, int bn, int m, int k, float* _scale_all ) {

	int N = bn + m - 1;
	//Eigen default to ColMajor: which means the array [N, k] is k rows and N columns
	auto rx_view = Map<const MatrixXcf> ((const complex<float>*)_rx_view, k, N);
	auto sxy_seg = Map<const MatrixXcf> ((const complex<float>*)_sxy_seg, k, bn);
	auto scale_all = Map<MatrixXcf> ((complex<float>*)_scale_all, k, bn);

	auto init_rx_blk = rx_view.leftCols(m);
	MatrixXcf sxxmat = (init_rx_blk * init_rx_blk.adjoint());
	
	auto sxx_inv = sxxmat.ldlt();
	for (int i = 0; i < bn; ++i) {
		scale_all.col(i) = sxx_inv.solve(sxy_seg.col(i).conjugate()).conjugate();
		if (i != bn - 1) {
			sxx_inv.rankUpdate(rx_view.col(i + m), 1);
			sxx_inv.rankUpdate(rx_view.col(i), -1);
		}
	}
}


struct ldlt_rls: Eigen::LDLT<MatrixXcd> {
	using super = Eigen::LDLT<MatrixXcd>;
	double lam;
	VectorXcd c;
	ldlt_rls (const ldlt_rls& other)
	: super(other.reconstructedMatrix()), lam(other.lam), c(other.c) { }

	//u_all layout (ch, step) (Eigen is colmajor)
	ldlt_rls (Map<MatrixXcd> & u_all, Map<VectorXcd> &d_all,
			double lam, const ldlt_rls* prev) : super(), lam(lam) {
		Index u_cols = u_all.cols(), u_rows = u_all.rows();
		MatrixXcd Ruu;
		VectorXcd Rud;
		double loglam = log(lam);
		double atten_prev = exp(loglam * u_cols);
		if (prev != nullptr) {
			assert(prev->ch() == u_rows);
			Ruu = prev->reconstructedMatrix() * atten_prev;
			Rud = Ruu * prev->c;
		} else {
			Ruu = MatrixXcd::Identity(u_rows, u_rows) * atten_prev;
			Rud.setZero(u_rows);
		}
		for (Index i = 0; i < u_cols; i++) {
			double atten = exp(loglam * (u_cols - i - 1));
			VectorXcd ua =  atten * u_all.col(i);
			Ruu +=  ua * u_all.col(i).adjoint();
			Rud += ua * conj(d_all(i));
		}
		super::compute(Ruu);
		c = solve(Rud);
	}
	void rls_get_c(Map<VectorXcd> & c1) const  { c1 = c; }
	void rls_update(Map<VectorXcd> &u, complex<double> err) {

		//update C += -ke*
		VectorXcd pu = solve(u);
		double upu = (u.adjoint() * pu).value().real();
		c += (conj(err) / -(lam + upu)) * pu;
		//update Ruu: attenuate by lam then rankUpdate
		for (Index i = 0, ii = ch(); i < ii; i++) {
			super::m_matrix.coeffRef(i, i) *= lam;
		}
		super::rankUpdate(u, 1);
	}
	

	Index ch() const { return super::m_matrix.cols(); }
	auto map_vector_ch(double *_u) const {
		return Map<VectorXcd> ((complex<double>*)_u, ch(), 1);
	}
};

//interface for python bindings of class LDLT_RLS

extern"C" void* new_ldlt_rls(
		int ch, int step, double *_u_all, double *_d_all,
		double lam, const void* prev) {
	Map<MatrixXcd> u_all((complex<double>*)_u_all, ch, step);
	Map<VectorXcd> d_all((complex<double>*)_d_all, step, 1);
	return new ldlt_rls(u_all, d_all, lam, (ldlt_rls*) prev);
}

extern "C" void free_ldlt_rls(void* ptr) {
	delete (ldlt_rls*) ptr;
}

extern "C" void* deepcopy_ldlt_rls(void* ptr) {
	return new ldlt_rls(*((ldlt_rls*)ptr));
}
extern "C" void ldlt_rls_get_c(void *ptr, double* _write) {
	auto& ldlt =  *((ldlt_rls*)ptr);
	auto c1 = ldlt.map_vector_ch(_write);
	ldlt.rls_get_c(c1);
}
extern "C" void ldlt_rls_update(void *ptr, double *_u, double* _err) {
	auto& ldlt =  *((ldlt_rls*)ptr);
	auto u = ldlt.map_vector_ch(_u);
	ldlt.rls_update(u, *(complex<double>*)_err);
}


