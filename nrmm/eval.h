#ifndef EVAL_H_
#define EVAL_H_

#include <algorithm>
#include <vector>
#include <map>

namespace npbayes
{
	double AdjustedRandIndex(const std::vector<int> &l0, const std::vector<int> &l1);
	void ComputeAvgLabels(const std::vector<std::vector<int>> labels, 
		std::vector<int> & avg_labels);
}
#endif