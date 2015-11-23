#include "eval.h"

namespace npbayes
{
	double AdjustedRandIndex(const std::vector<int> &l0, const std::vector<int> &l1)
	{
		int n = l0.size();
		int k0 = *std::max_element(l0.begin(), l0.end()) + 1;
		int k1 = *std::max_element(l1.begin(), l1.end()) + 1;
		std::vector<std::vector<int>> Cont(k0, std::vector<int>(k1, 0));
		for (int i = 0; i < n; ++i)
			Cont[l0[i]][l1[i]]++;

		std::vector<int> row_sum(Cont.size(), 0);
		std::vector<int> col_sum(Cont[0].size(), 0);
		for (int i = 0; i < Cont.size(); ++i) {
			for (int j = 0; j < Cont[i].size(); ++j) {
				row_sum[i] += Cont[i][j];
				col_sum[j] += Cont[i][j];
			}
		}

		double a = 0;
		for (int i = 0; i < k0; ++i) {
			for (int j = 0; j < k1; ++j) {
				if (Cont[i][j]) a += 0.5 * Cont[i][j] * (Cont[i][j] - 1);
			}
		}

		double b1 = 0;
		double b2 = 0;
		for (int i = 0; i < k0; ++i) {
			if (row_sum[i]) b1 += 0.5 * row_sum[i] * (row_sum[i] - 1);
		}
		for (int i = 0; i < k1; ++i) {
			if (col_sum[i]) b2 += 0.5 * col_sum[i] * (col_sum[i] - 1);
		}

		double c = 0.5 * n * (n - 1);

		return (a - b1*b2 / c) / (0.5*(b1 + b2) - b1*b2 / c);
	}

	void ComputeAvgLabels(const std::vector<std::vector<int>> labels, 
		std::vector<int> & avg_labels)
	{
		int m = labels.size();
		if (m == 0) return;
		int n = labels[0].size();
		std::map<std::pair<int, int>, int> cnt;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n - 1; ++j) {
				for (int k = j + 1; k < n; ++k) {
					if (labels[i][j] == labels[i][k]) 
						cnt[std::pair<int, int>(j, k)]++;					
				}
			}
		}
		avg_labels.assign(n, 0);
		std::vector<bool> unchecked(n, true);
		int thres = m / 2;
		int l = 0;
		for (int i = 0; i < n; ++i) {
			if (unchecked[i]) {
				unchecked[i] = false;
				avg_labels[i] = l;
				for (int j = i + 1; j < n; ++j) {
					if (cnt[std::pair<int, int>(i, j)] > thres) {
						unchecked[j] = false;
						avg_labels[j] = l;
					}
				}
				++l;
			}
		}
	}
}