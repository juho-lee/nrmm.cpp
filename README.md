# nrmm.cpp
MCMC inference algorithms for normalized random measure mixtures 
- Marginal Gibbs Sampler
- Split-Merge Sampler
- Tree-guided MCMC (NIPS 2015, Paper: http://arxiv.org/abs/1511.05650)

Coded & tested with Microsoft Visual Studio 2013 on Windows machine

Usage: in Relase folder, type
nrmm data output bm nrm init sampler params
- data: data name (e.g., toy, 10k, nips)
- output: output folder name 
- bm: base measure, 0 for NormalWishart and 1 for Multinomial-Dirichlet
- nrm: NRM, 0 for DP and 1 for NGGP
- init: initialization option, 0 for exact IBHC and 1 for noisy IBHC (see experimental section of the paper)
- sampler: sampler, 0 for Gibbs, 1 for split-merge and 2 for TGMCMC
- params: parameters for each samplers
  * nrmm data output bm nrm init 0 subset et_thres
  * nrmm data output bm nrm init 1 subset et_thres
    * subset for subset size (see paper) and et_thres for total running time
  * nrmm data output bm nrm init 2 num_sm depth et_thres
    * num_sm for the parameter G and depth for the parameter D in the paper
 
Run sample_script.bat and display_results.m to produce results for toy dataset

The nips data was accquired from https://archive.ics.uci.edu/ml/datasets/Bag+of+Words.



