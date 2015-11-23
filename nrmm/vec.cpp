#include "vec.h"

namespace npbayes
{
	int ReadMat(const char *filename, Mat &X)
	{
		std::vector< std::vector<double> > buff;
		std::string line;
		std::ifstream ifs(filename);
		if (!ifs) return 0;
		while (getline(ifs, line)) {
			std::istringstream is(line);
			buff.push_back(std::vector<double>(std::istream_iterator<double>(is),
				std::istream_iterator<double>()));
		}
		ifs.close();
		X.resize(buff.size(), buff[0].size());
		for (int i = 0; i < buff.size(); ++i) {
			for (int j = 0; j < buff[0].size(); ++j)
				X(i, j) = buff[i][j];
		}
		return 1;
	}

	int WriteMat(const char *filename, const Mat &X)
	{
		std::ofstream ofs(filename);
		if (!ofs) return 0;
		for (int i = 0; i < X.rows(); ++i) {
			for (int j = 0; j < X.cols() - 1; ++j)
				ofs << X(i, j) << " ";
			ofs << X(i, X.cols() - 1) << std::endl;
		}
		ofs.close();
		return 1;
	}
}