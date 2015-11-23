function [max_lj, et, ess] = compute_stats(data, sampler, ntrials, noess)
    if nargin == 3 
        noess = 0;
    end
    path = sprintf('%s/%s', data, sampler);    
    max_lj = zeros(ntrials, 1);
    et = zeros(ntrials, 1);
    ess = zeros(ntrials, 1);
    for i = 1 : ntrials
        S = dlmread(sprintf('%s_%d/log.txt', path, i-1));
        max_lj(i) = max(S(:, 3));
        et(i) = mean(S(2:end,1) - S(1:end-1,1));
        if ~noess
            ess(i) = dlmread(sprintf('%s_%d/ess.txt', path, i-1));
        end
    end
    
    if ~noess
        fprintf('%.6f (%.6f), %.6f (%.6f), %.6f (%.6f)\n', mean(max_lj), std(max_lj), ...
        mean(et), std(et), mean(ess), std(ess));
    else
        fprintf('%.6f (%.6f), %.6f (%.6f)\n', mean(max_lj), std(max_lj), mean(et), std(et));
    end