function X = gen_mixture(gt_labels, dist, varargin)    
    n = length(gt_labels);
    k = max(gt_labels);
    switch dist
        case 'multinomial'
            d = 20;
            m = 10;
            a = 0.5;
            while ~isempty(varargin)
                switch lower(varargin{1})  
                    case 'd'
                        d = varargin{2};
                        varargin(1:2) = [];
                    case 'm'
                        m = varargin{2};
                        varargin(1:2) = [];
                    case 'a'
                        a = varargin{2};
                        varargin(1:2) = [];
                    otherwise
                        fprintf('Wrong option: %s\n', varargin{1});
                        varargin(1:2) = [];
                end
            end
            X = zeros(d, n);
            for j = 1 : k
                w = rand_dir(a*ones(1,d));
                X(:,gt_labels==j) = mnrnd(m, w, sum(gt_labels==j))';
            end            
        case 'gaussian'
            d = 2;
            r = 0.08;
            df = 2*d;
            while ~isempty(varargin)
                switch lower(varargin{1})  
                    case 'd'
                        d = varargin{2};
                        varargin(1:2) = [];
                    case 'r'
                        r = varargin{2};
                        varargin(1:2) = [];
                    case 'df'
                        df = varargin{2};
                        varargin(1:2) = [];
                    otherwise
                        fprintf('Wrong option: %s\n', varargin{1});
                        varargin(1:2) = [];
                end
            end                 
            X = zeros(d, n);
            Scale = rand(d);
            Scale = Scale*Scale' + eye(d);
            cholScale = chol(Scale)';            
            for j = 1 : k
                [~, cholLambda] = rand_wishart(df,[],cholScale);
                cholSigma = (cholLambda^-1)';
                mu = cholSigma*randn(d,1)/sqrt(r);
                nj = sum(gt_labels==j);
                X(:,gt_labels==j) = repmat(mu,1,nj) + cholSigma*randn(d,nj);
            end
        otherwise
            fprintf('Unsupported distribution: %s\n', dist);
            X = [];
    end