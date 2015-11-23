function [labels, k] = gen_dp(n, alpha)
    labels = zeros(n, 1);
    k = 0;
    cnt = [];
    for i = 1 : n
        w = [cnt alpha];
        j = rand_categorial(w);
        labels(i) = j;
        if (j <= k)
            cnt(j) = cnt(j) + 1;
        else
            cnt = cat(2, cnt, k+1);
            k = k + 1;
        end
    end
