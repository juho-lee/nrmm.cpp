#include "normal_wishart.h"

namespace npbayes
{
	void NormalWishartSuffStats::Add(SuffStats *ss)
	{
		NormalWishartSuffStats *nwsss = (NormalWishartSuffStats *)(ss);
		n_ += nwsss->n_; 
		x_ += nwsss->x_;
		X_ += nwsss->X_;
	}

	void NormalWishartSuffStats::Subtract(SuffStats *ss)
	{
		NormalWishartSuffStats *nwsss = (NormalWishartSuffStats *)(ss);
		n_ -= nwsss->n_;
		x_ -= nwsss->x_;
		X_ -= nwsss->X_;
	}

	NormalWishart::NormalWishart(const Mat &X, double c)
	{
		xi_ = X.rowwise().mean();
		int n = X.cols();
		Mat Psi = X*X.transpose() / (n - 1) - n*xi_*xi_.transpose() / (n - 1);
		d_ = xi_.rows();
		Psi /= pow(Psi.determinant() / c, 1.0 / d_);
		r_ = 0.1;
		nu_ = d_ + 6;
		chol_Psi_inv_ = Psi.inverse().llt().matrixL();
		Psi_rxixit_ = Psi + (r_*xi_)*xi_.transpose();
		nu_log_det_Psi_ = 0.5*nu_*log(Psi.determinant());
	}

	double NormalWishart::LogMvnGammaRatio(int nu, int n, int d)
	{
		if (n == 1 && d == 2)
			return log(0.5*(nu - 1));
		int half_n = n / 2; 
		double val = -0.69314718055994529*d*half_n;
		if (n % 2 == 0) {
			for (int i = 1; i <= d; i++) {
				for (int j = 1; j <= half_n; j++)
					val += log((double)(nu - 1 - i + 2 * j));
			}
		}
		else {
			for (int i = 1; i <= d; i++) {
				for (int j = 1; j <= half_n; j++)
					val += log((double)(nu - i + 2 * j));
				val += LogGamma(0.5 * (double)(nu + 2 - i))
					- LogGamma(0.5 * (double)(nu + 1 - i));
			}
		}
		return val;
	}

	double NormalWishart::LogMarginal(SuffStats *ss) const
	{
		const static double log_PI = log(PI);
		NormalWishartSuffStats *nwss = (NormalWishartSuffStats*)(ss);		
		double r_n = r_ + nwss->n_;
		int nu_n = nu_ + nwss->n_;
		Vec xi_n = (r_*xi_ + nwss->x_) / r_n;
		Mat Psi_n = Psi_rxixit_ + nwss->X_ - r_n*xi_n*xi_n.transpose();
		return (-0.5*nwss->n_*d_*log_PI
			+ 0.5*d_*(log(r_) - log(r_n))
			+ nu_log_det_Psi_ - 0.5*nu_n*log(Psi_n.determinant())
			+ LogMvnGammaRatio(nu_, nwss->n_, d_));
	}

	double NormalWishart::LogMarginal(SuffStats *ss0, SuffStats *ss1) const
	{
		const static double log_PI = log(PI);
		NormalWishartSuffStats *nwss0 = (NormalWishartSuffStats*)(ss0);
		NormalWishartSuffStats *nwss1 = (NormalWishartSuffStats*)(ss1);
		int n = nwss0->n_ + nwss1->n_;
		double r_n = r_ + n;
		int nu_n = nu_ + n;
		Vec xi_n = (r_*xi_ + nwss0->x_ + nwss1->x_) / r_n;
		Mat Psi_n = Psi_rxixit_ + nwss0->X_ + nwss1->X_ - r_n*xi_n*xi_n.transpose();
		return (-0.5*n*d_*log_PI
			+ 0.5*d_*(log(r_) - log(r_n))
			+ nu_log_det_Psi_ - 0.5*nu_n*log(Psi_n.determinant())
			+ LogMvnGammaRatio(nu_, n, d_));
	}
}