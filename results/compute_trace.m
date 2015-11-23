function [lj_mean, lj_std, nc_mean, nc_std, et] ...
    = compute_trace(data, sampler, ntrials)
    et_cell = cell(1, ntrials);
    lj_cell = cell(1, ntrials);
    nc_cell = cell(1, ntrials);
    path = sprintf('%s/%s', data, sampler);
    
    timestep = [];    
    max_time = -Inf;
    for i = 1 : ntrials
        S = dlmread(sprintf('%s_%d/log.txt', path, i-1));
        et_cell{i} = S(:, 1);
        nc_cell{i} = S(:, 2);
        lj_cell{i} = S(:, 3);
        timestep = mean([timestep mean(et_cell{i}(2:end)-et_cell{i}(1:end-1))]);        
        max_time = max(max_time, et_cell{i}(end));
    end
        
    et = 0 : timestep : max_time;
    lj_mean = zeros(1, length(et));
    lj_std = zeros(1, length(et));
    nc_mean = lj_mean;
    nc_std = lj_mean;
    for i = 1 : length(et)    
        lj = [];        
        nc = [];
        for j = 1 : ntrials
            ind = find(et_cell{j}>et(i), 1);
            lj = cat(1, lj, lj_cell{j}(ind));  
            nc = cat(1, nc, nc_cell{j}(ind));
        end        
        lj_mean(i) = mean(lj);       
        lj_std(i) = std(lj);
        nc_mean(i) = mean(nc);
        nc_std(i) = mean(nc);
    end
