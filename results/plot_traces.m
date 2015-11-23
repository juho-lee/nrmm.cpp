function plot_traces(data, sampler, label, ntrials, maxtime, bias)    
    if nargin == 5
        bias = 0.01;
    end
    colors = gen_colors(length(sampler));
    hold on;
    maxval = -inf;
    for i = 1 : length(sampler)
        [lj, ~, ~, ~, et] = compute_trace(data, sampler{i}, ntrials);
        maxval = max(max(lj), maxval);
        plot(et, lj, 'linewidth', 2, 'color', colors(i,:));
    end
    set(gca, 'fontsize', 20);
    xlim([0, maxtime]);
    ylim([-inf, maxval + abs(bias*maxval)]);
    xlabel('time [sec]', 'fontsize', 30);
    ylabel('average log-likelihood', 'fontsize', 30);
    legend(label, 'fontsize', 30, 'location', 'southeast');
end