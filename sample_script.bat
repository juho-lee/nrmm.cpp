%for /l %%x in (0, 1, 4) do (nrmm toy nggp_gibbs_%%x 0 1 1 0 0 10)
%for /l %%x in (0, 1, 4) do (nrmm toy nggp_sm_%%x 0 1 1 1 0 10)
%for /l %%x in (0, 1, 4) do (nrmm toy nggp_tgmcmc_%%x 0 1 1 2 10 2 10)

for /l %%x in (0, 1, 4) do (nrmm toy dp_gibbs_%%x 0 0 1 0 0 10)
for /l %%x in (0, 1, 4) do (nrmm toy dp_sm_%%x 0 0 1 1 0 10)
for /l %%x in (0, 1, 4) do (nrmm toy dp_tgmcmc_%%x 0 0 1 2 10 2 10)

nrmm 10k dp_gibbs 0 0 1 0 0 200
nrmm 10k dp_tgmcmc 0 0 1 2 10 2 200
nrmm 10k nggp_gibbs 0 1 1 0 0 200
nrmm 10k nggp_tgmcmc 0 1 1 2 10 2 200