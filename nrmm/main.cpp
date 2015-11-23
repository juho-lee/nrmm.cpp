#include "dp.h"
#include "nggp.h"
#include "normal_wishart.h"
#include "mult_dir.h"
#include "nrmm_gibbs_sampler.h"
#include "nrmm_split_merge_sampler.h"
#include "tgmcmc.h"

using namespace npbayes;

int InitNormalWishart(const char *filename,
	std::vector<SuffStats*> &ss, BaseMeasure *&H)
{
	Mat X;
	if (ReadMat(filename, X)) {
		X.transposeInPlace();
		H = new NormalWishart(X);
		ss.assign(X.cols(), 0);
		for (int i = 0; i < ss.size(); ++i)
			ss[i] = new NormalWishartSuffStats(X.col(i));
		return 1;
	}
	else return 0;
}

int InitMultDir(const char *filename, double alpha,
	std::vector<SuffStats*> &ss, BaseMeasure *&H)
{
	std::ifstream ifs(filename);
	if (!ifs) return 0;
	int num_docs, num_words, dump;
	ifs >> num_docs >> num_words >> dump;
	std::vector<std::map<int, int>> hist(num_docs, std::map<int, int>());
	int docid, wordid, cnt;
	while (ifs >> docid >> wordid >> cnt)
		hist[docid - 1][wordid - 1] += cnt;
	ss.assign(num_docs, 0);
	for (int i = 0; i < num_docs; ++i)
		ss[i] = new MultDirSuffStats(hist[i]);
	H = new MultDir(num_words, alpha);
	return 1;
}

int main(int argc, char **argv)
{
#ifdef _DEBUG
	char *data = "nips";
	char *output = "temp";
	int bmopt = 1;
	int nrmopt = 1;
	int initopt = 0;
	int samopt = 2;
	double et_thres = 1000;
#else
	if (argc != 10 && argc != 9) {
		std::cout << "nrmm data output bmopt nrmopt initopt samopt params" << std::endl;
		std::cout << "bmopt: 0 (NormalWishart), 1 (MultDir)" << std::endl;
		std::cout << "nrmopt: 0 (DP), 1 (NGGP)" << std::endl;
		std::cout << "initopt: 0 (Exact IBHC), 1 (Noisy IBHC)" << std::endl;
		std::cout << "samopt: 0 (Gibbs), 1 (Split-Merge), 2 (TGMCMC)" << std::endl;
		std::cout << "nrmm data output bmopt nrmopt initopt 0 subset et_thres " << std::endl;
		std::cout << "nrmm data output bmopt nrmopt initopt 1 subset et_thres" << std::endl;
		std::cout << "nrmm data output bmopt nrmopt initopt 2 num_sm depth et_thres" << std::endl;
		return 0;
	}
	srand(time(0));
	char *data = argv[1];
	char *output = argv[2];
	int bmopt = atoi(argv[3]);
	int nrmopt = atoi(argv[4]);
	int initopt = atoi(argv[5]);
	int samopt = atoi(argv[6]);
	double et_thres;
	if (samopt == 0 || samopt == 1) et_thres = atof(argv[8]);
	else et_thres = atof(argv[9]);
#endif
	char path[BUFFSIZE];
	sprintf_s<BUFFSIZE>(path, "%s/%s/%s", resultspath, data, output);
	char buff[BUFFSIZE];
	sprintf_s<BUFFSIZE>(buff, "mkdir \"%s\"", path);
	system(buff);

	sprintf_s<BUFFSIZE>(buff, "%s/settings.txt", path);
	std::ofstream ofs(buff);
	ofs << "data: " << data << std::endl;

	std::vector<SuffStats*> ss;
	BaseMeasure *H = 0;
	if (bmopt == 0) {
		std::cout << "Base Measure: Normal Wishart" << std::endl;
		ofs << "Base Measure: Normal Wishart" << std::endl;
		sprintf_s<BUFFSIZE>(buff, "%s/normal_wishart/%s/X.txt", datapath, data);
		if (!InitNormalWishart(buff, ss, H)) {
			std::cout << "Error: data not available." << std::endl;
			ofs.close();
			return 0;
		}
	}
	else if (bmopt == 1) {
		std::cout << "Base Measure: Mult Dir" << std::endl;
		ofs << "Base Measure: Mult Dir" << std::endl;
		sprintf_s<BUFFSIZE>(buff, "%s/mult_dir/%s/docword.txt", datapath, data);
		double alpha = 0.1;
		if (!InitMultDir(buff, alpha, ss, H)) {
			std::cout << "Error: data not available." << std::endl;
			ofs.close();
			return 0;
		}
	}
	else {
		std::cout << "Error: unsupported base measure." << std::endl;
		ofs.close();
		return 0;
	}

	NRM *mu = 0;
	if (nrmopt == 0) {
		std::cout << "NRM: DP" << std::endl;
		ofs << "NRM: DP" << std::endl;
		mu = new DP(H);
	}
	else if (nrmopt == 1) {
		std::cout << "NRM: NGGP" << std::endl;
		ofs << "NRM: NGGP" << std::endl;
		mu = new NGGP(H);
	}
	else {
		std::cout << "Error: unsupported NRM." << std::endl;
		ofs.close();
		return 0;
	}

	NRMMSampler *sampler = 0;
	if (samopt == 0) {
		std::cout << "Sampler: Gibbs" << std::endl;
		ofs << "Sampler: Gibbs" << std::endl;
		ofs << "Subset Size: " << atoi(argv[7]) << std::endl;
		sampler = new NRMMGibbsSampler(mu, ss, atoi(argv[7]));
	}
	else if (samopt == 1) {
		std::cout << "Sampler: Split-Merge" << std::endl;
		ofs << "Sampler: Split-Merge" << std::endl;
		ofs << "Subset Size: " << atoi(argv[7]) << std::endl;
		sampler = new NRMMSplitMergeSampler(mu, ss, atoi(argv[7]));
	}
	else if (samopt == 2) {
		std::cout << "Sampler: TGMCMC" << std::endl;
		ofs << "Sampler: TGMCMC" << std::endl;

		ofs << "Number of SM: " << atoi(argv[7]) << std::endl;
		ofs << "Depth: " << atoi(argv[8]) << std::endl;
		sampler = new TGMCMC(mu, ss, atoi(argv[7]), atoi(argv[8]));
	}
	else {
		std::cout << "Error: unsupported sampler." << std::endl;
		ofs.close();
		return 0;
	}
	ofs.close();

	sprintf_s<BUFFSIZE>(buff, "%s/log.txt", path);
	ofs.open(buff);
	ofs.precision(PRECISION);
	std::cout.precision(PRECISION);
	clock_t bc = clock();
	if (samopt == 0 || samopt == 1) {
		TGMCMC *temp = new TGMCMC(mu, ss);
		temp->Init(bool(initopt));
		std::vector<int> labels;
		temp->GetLabels(labels);
		sampler->Init(labels);
	}
	else sampler->Init(bool(initopt));
	double et_cum = ElapsedSecs(bc);
	double lj = sampler->LogJoint();
	int nc = sampler->NumClusters();

	std::cout << "initialization: ";
	if (initopt) std::cout << "noisy IBHC" << std::endl;
	else std::cout << "exact IBHC" << std::endl;

	std::cout << "init: (" << et_cum << " secs), "
		<< nc << " clusters, " << lj << std::endl;
	ofs << et_cum << " " << nc << " " << lj << std::endl;

	int iter = 1;
	while (et_cum < et_thres) {	
		bc = clock();
		sampler->Sweep();
		et_cum += ElapsedSecs(bc);
		lj = sampler->LogJoint();
		nc = sampler->NumClusters();
		std::cout << "iter " << iter++ << " (" << et_cum << " secs), "		
			<< nc << " clusters, " << lj << std::endl;
		ofs << et_cum << " " << nc << " " << lj << std::endl;
	}
	ofs.close();

	sprintf_s<BUFFSIZE>(buff, "%s/addinfo.txt", path);
	ofs.open(buff);
	sampler->PrintAdditionalInfo(ofs);
	ofs.close();

	return 0;
}
