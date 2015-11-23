%%
clc;
clear;

%%
data = 'toy';
ntrials = 5;
sampler = {'dp_gibbs';'dp_sm';'dp_tgmcmc'};
leg = {'Gibbs';'SM';'TGMCMC'};

colors = ['r';'g';'b'];
maxtime = 10;

%%
hold on;
maxval = -inf;
for i = 1 : length(sampler)
    [lj, ~, ~, ~, et] = compute_trace(data, sampler{i}, ntrials);
    plot(et, lj, 'color', colors(i,:), 'linewidth', 1.5);
    maxval = max(maxval, max(lj(:)));
end
legend(leg, 'location', 'southeast', 'fontsize', 30);
xlim([0 maxtime]);
ylim([-inf maxval]);
set(gca, 'fontsize', 20);
xlabel('time [sec]', 'fontsize', 20);
ylabel('log-likelihood', 'fontsize', 20);
