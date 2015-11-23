#include "nggp.h"

namespace npbayes
{
	const double NGGP::a_alpha_ = 1.0;
	const double NGGP::b_alpha_ = 1.0;
	const double NGGP::a_sigma_ = 1.5;
	const double NGGP::b_sigma_ = 1.5;

	double NGGP::LogKappa(int n) const
	{
		return log(alpha_) + log(sigma_) + LogGamma(n - sigma_)
			- (n - sigma_)*log(1 + u_) - LogGamma(1 - sigma_);
	}

	double NGGP::LogKappaJoin(int n) const
	{
		return log(n - sigma_);
	}

	double NGGP::LogKappaNew(void) const
	{
		return log(alpha_) + log(sigma_) + sigma_*log(1 + u_);
	}

	double NGGP::LogJointConst(const ClusterSet &clset) const
	{
		int n = 0;
		foreach(it, clset) n += it->n;
		return (n - 1)*log(u_) - alpha_*(pow(1 + u_, sigma_) - 1) - LogGamma(n);
	}

	double LogFuGibbs(double x, double alpha, double sigma, int n, int k)
	{
		double ex = exp(x);
		return x*n - alpha*pow(1 + ex, sigma) - (n - k*sigma)*log(1 + ex);
	}

	void NGGP::SampleU(const ClusterSet &clset)
	{
		int n = 0;
		foreach(it, clset) n += it->n;
		fun log_f = std::bind(LogFuGibbs, _1, alpha_, sigma_, n, clset.size());
		u_ = exp(SliceSample(log_f, log(u_), 0.001));
	}

	double LogFsigmaGibbs(double x, double u, double alpha,
		const ClusterSet &clset)
	{
		int k = clset.size();
		double y = (NGGP::a_sigma_ + k - 1)*log(x)
			+ (NGGP::b_sigma_ - 1)*log(1 - x) - alpha*pow(1 + u, x);
		foreach(it, clset) {
			y += LogGamma(it->n - x) - LogGamma(1 - x)
				- (it->n - x) * log(1 + u);
		}
		return y;
	}

	void NGGP::SampleHyperparams(const ClusterSet &clset)
	{
		int n = 0, k = clset.size();
		foreach(it, clset) n += it->n;
		alpha_ = RandGamma(a_alpha_ + k, b_alpha_ + pow(1 + u_, sigma_) - 1);
		fun log_f = std::bind(LogFsigmaGibbs, _1, u_, alpha_, clset);
		sigma_ = SliceSample(log_f, sigma_, 0.01, 0.99, 0.01);
	}

	//double LogFuSlice(double x, int n,
	//	double alpha, double sigma, double w_sum, double L)
	//{
	//	double ex = exp(x);
	//	return (x*n - ex*w_sum - alpha*pow(1 + ex, sigma)
	//		*(1 - Q(-sigma, (1 + ex)*L)));
	//}

	//void NGGP::SampleU(double L, const StickSet &stkset)
	//{
	//	int n = 0;
	//	double w_sum = 0;
	//	foreach(it, stkset) {
	//		w_sum += it->w;
	//		n += it->n;
	//	}
	//	fun log_f = std::bind(LogFuSlice, _1, n, alpha_, sigma_, w_sum, L);
	//	u_ = exp(SliceSample(log_f, log(u_), 0.001));
	//}

	//double LogFsigmaSlice(double x, double u, double alpha,
	//	int k, double L, double log_w_sum)
	//{
	//	return ((NGGP::a_sigma_ + k - 1) * log(x) - k*LogGamma(1 - x)
	//		+ (NGGP::b_sigma_ - 1)*log(1 - x) - x*log_w_sum
	//		- alpha*pow(1 + u, x)*(1 - Q(-x, (1 + u)*L)));
	//}

	//void NGGP::SampleHyperParams(double L, const StickSet &stkset)
	//{
	//	int k = stkset.size();
	//	alpha_ = RandGamma(a_alpha_ + k, b_alpha_ + pow(1 + u_, sigma_)
	//		*(1 - Q(-sigma_, (1 + u_)*L)) - 1);
	//	double log_w_sum = 0;
	//	foreach(it, stkset) log_w_sum += log(it->w);
	//	fun log_f = std::bind(LogFsigmaSlice, _1, u_, alpha_, k, L, log_w_sum);
	//	sigma_ = SliceSample(log_f, sigma_, 0.01, 0.99, 0.01);
	//}

	//void NGGP::SampleJumps(StickSet &stkset)
	//{
	//	for (auto it = stkset.begin(); it != stkset.end();) {
	//		if ((*it)->n == 0) {
	//			delete *it;
	//			it = stkset.erase(it);
	//		}
	//		else {
	//			(*it)->w = RandGamma((*it)->n - sigma_, 1 + u_);				
	//			++it;
	//		}
	//	}
	//}

	//void NGGP::SampleJumps(double L, StickSet &stkset, int maxlen)
	//{
	//	double t = L;
	//	int len = 0;
	//	while (true) {
	//		double log_r = log(RandExp(1.0));
	//		double log_thres = log(alpha_) + log(sigma_)
	//			- (1 + sigma_)*log(t) - (1 + u_)*t
	//			- log(1 + u_) - LogGamma(1 - sigma_);
	//		if (log_r > log_thres) break;
	//		else {
	//			double tp = t - log(1 - exp(log_r - log_thres)) / (1 + u_);
	//			if (Randu() < exp(-(1 + sigma_)*(log(tp) - log(t)))) {
	//				stkset.insert(new Stick(tp));
	//				if (++len > maxlen) break;
	//			}
	//			t = tp;				
	//		}
	//	}
	//}
}