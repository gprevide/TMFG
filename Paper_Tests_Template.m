addpath('C:\Users\Guido\Documents\Bauhaus\PhD\matlab_bgl', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\pmfg', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\BCT', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\RANDRAW')
    



%%%
% POWER DISTRIBUTION TEST
%%%
for alpha = [1 2]
    %----------------------|
    fprintf('Pareto Alpha %d\n', alpha)
    fprintf('N|PMFG|TMFG|Ratio|PlanarEdges|Totaledges|Edgeratio|Avg|Stdev|Skew\n');
    for n = [ 500 700 1000 1200 1500]
        x = randraw('pareto', [1, alpha], [1 n*(n-1)/2]);
        % x = randraw('exp', lambda, [1 n*(n-1)/2]);
        % x = randraw('beta', [alpha 2], [1 n*(n-1)/2]);
        x = (x - min(x)) ./ (max(x) - min(x));
        W = triu(ones(n,n),1);
        W(W==1) = x;
        r = W .* W;
        r = max(r,r');
        fprintf('%d|', n);
        tic;[TT, cliques, triangles, tc  ] = TMFG_mc(r); elapsed = toc;
        fprintf('TMFG: %f, %f|', sum(TT(:))+0.0, elapsed);
        tic;[TT1, cliques, triangles, tc  ] = TMFGT1_mc(r); elapsed = toc;
        fprintf('TMFGT1: %f, %f|', sum(TT1(:))+0.0, elapsed);
        tic;PP = pmfg(r);elapsed = toc;
        fprintf('PMFG: %f, %f|', sum(PP(:))+0.0, elapsed);
%         tic;[TT21, cliques, separators, peo ] = TMFG_v2_1(r); elapsed = toc;
%         fprintf('TMFG_V2_1: %f, %f\n', sum(TT21(:))+0.0, elapsed);
        tic;[TT22, cliques, separators, peo ] = TMFG_v2_2(r); elapsed = toc;
        fprintf('TMFG_V2_2: %f, %f|', sum(TT22(:))+0.0, elapsed);
        tic;[TT23, cliques, separators, peo ] = TMFG_v2_3(r); elapsed = toc;
        fprintf('TMFG_V2_3: %f, %f\n', sum(TT23(:))+0.0, elapsed);
    end
end

%%%
% BETA DISTRIBUTION TEST
%%%
% for params = [  2 3 4 5 10;   4 3 2 1 0.5]
%     fprintf('Beta distribution alpha: %d beta: %d\n', params(1), params(2))
%     for n = [ 100 200 300 400 500 ]
%         % W = rand(n,n); % Uniform
%         % x = randraw('pareto', [1, alpha], [1 n*(n-1)/2]);
%         % x = randraw('exp', lambda, [1 n*(n-1)/2]);
%         x = randraw('beta', [params(1) params(2)], [1 n*(n-1)/2]);
%         x = x ./ max(x);
%         W = triu(ones(n,n),1);
%         W(W==1) = x;
%         r = W .* W;
%         r = max(r,r');
%         fprintf('%d|', n);
%         tic;[TT, cliques, triangles, tc  ] = TMFG_mc(r); elapsed = toc;
%         fprintf('TMFG: %f, %f|', sum(TT(:))+0.0, elapsed);
%         tic;[TT1, cliques, triangles, tc  ] = TMFGT1_mc(r); elapsed = toc;
%         fprintf('TMFGT1: %f, %f|', sum(TT1(:))+0.0, elapsed);
%         tic;PP = pmfg(r);elapsed = toc;
%         fprintf('PMFG: %f, %f|', sum(PP(:))+0.0, elapsed);
%         % tic;[TT21, cliques, separators, peo ] = TMFG_v2_1(r); elapsed = toc;
%         % fprintf('TMFG_V2_1: %f, time: %f\n', sum(TT21(:))+0.0, elapsed);
%         tic;[TT22, cliques, separators, peo ] = TMFG_v2_2(r); elapsed = toc;
%         fprintf('TMFG_V2_2: %f, %f|', sum(TT22(:))+0.0, elapsed);
%         tic;[TT23, cliques, separators, peo ] = TMFG_v2_3(r); elapsed = toc;
%         fprintf('TMFG_V2_3: %f, %f\n', sum(TT23(:))+0.0, elapsed);
%     end
% end

%%%
% LOGNORMAL DISTRIBUTION TEST
%%%
% fprintf('Lognormal\n')
% fprintf('N|PMFG|PTIME|TMFG|TTIME|Ratio|PlanarEdges|Totaledges|Edgeratio|Avg|Stdev|Skew\n');
% for n = [ 100 200 300 400 500 700 900 1000 1200 1500]
%     % W = rand(n,n); % Uniform
%     % x = randraw('pareto', [1, alpha], [1 n*(n-1)/2]);
%     % x = randraw('exp', lambda, [1 n*(n-1)/2]);
%     x = randraw('lognorm', [], [1 n*(n-1)/2]);
%     x = x ./ max(x);
%     W = triu(ones(n,n),1);
%     W(W==1) = x;
%     W = max(W, W');
%     tic; p = pmfg(W); pte = toc;
%     % square to get TMFG weights
%     W = W .* W;
%     tic; [t c s peo] = TMFG(W); tte = toc;
%     t = sqrt(t);
%     fprintf('%04d|%9.2f|%f|%9.2f|%f|%f|%d|%d|%f|', n, sum(p(:))+0, pte, sum(t(:))+0, tte, (sum(t(:))+0)/(sum(p(:))+0), 3*(n-2) , n*(n-1)/2, (3*(n-2))/(n*(n-1)/2) );
%     fprintf('%f|%f|%f\n', mean(x), std(x), skewness(x));
% end

%%%
% UNIFORM DISTRIBUTION TEST
%%%
fprintf('Uniform\n')
fprintf('Time test\n');
test = 0;
for n = [4000 8000 16000 32000]
    fprintf('%d|', n);
    x = randraw('uniform', [0 1], [1 n*(n-1)/2]);
    x = x ./ max(x);
    W = triu(ones(n,n),1);
    W(W==1) = x;
    W = max(W, W');
    r = W .* W;
    test = test + 1;
    tic;[TT, cliques, triangles, tc  ] = TMFG_mc(r); elapsed = toc;
    times(test, 1) = elapsed;
    fprintf('TMFG| %f| %f|', sum(TT(:))+0.0, elapsed);
    tic;[TT1, cliques, triangles, tc  ] = TMFGT1_mc(r); elapsed = toc;
    fprintf('TMFGT1| %f| %f\n', sum(TT1(:))+0.0, elapsed);
    times(test, 2) = elapsed;
%     tic;PP = pmfg(r);elapsed = toc;
%     fprintf('PMFG| %f| %f\n', sum(PP(:))+0.0, elapsed);
%     times(test, 3) = elapsed;
    % tic;[TT21, cliques, separators, peo ] = TMFG_v2_1(r); elapsed = toc;
    % fprintf('TMFG_V2_1: %f, time: %f\n', sum(TT21(:))+0.0, elapsed);
%     tic;[TT22, cliques, separators, peo ] = TMFG_v2_2(r); elapsed = toc;
%     fprintf('TMFG_V2_2| %f| %f\n', sum(TT22(:))+0.0, elapsed);
%     tic;[TT23, cliques, separators, peo ] = TMFG_v2_3(r); elapsed = toc;
%     fprintf('TMFG_V2_3: %f, %f\n', sum(TT23(:))+0.0, elapsed);
end
times;

%%%
% RANDOM CORRELATION MATRICES TEST (FACTORS METHOD)
%%%
% fprintf('Random corr matrices (Factors method)\n')
% fprintf('N|PMFG|PTIME|TMFG|TTIME|Ratio|PlanarEdges|Totaledges|Edgeratio|Avg|Stdev|Skew\n');
% for n = [50 100 150 200 300 500 750 1000]
%     W = gencorr(n, 'factors'); W = (W + W')/2; % get rid of minor differences
%     r = W .* W;
%     r = max(r,r');
%     fprintf('%d|', n);
%     tic;[TT, cliques, triangles, tc  ] = TMFG_mc(r); elapsed = toc;
%     fprintf('TMFG: %f, %f|', sum(TT(:))+0.0, elapsed);
%     tic;[TT1, cliques, triangles, tc  ] = TMFGT1_mc(r); elapsed = toc;
%     fprintf('TMFGT1: %f, %f|', sum(TT1(:))+0.0, elapsed);
%     tic;PP = pmfg(r);elapsed = toc;
%     fprintf('PMFG: %f, %f|', sum(PP(:))+0.0, elapsed);
%     % tic;[TT21, cliques, separators, peo ] = TMFG_v2_1(r); elapsed = toc;
%     % fprintf('TMFG_V2_1: %f, time: %f\n', sum(TT21(:))+0.0, elapsed);
%     tic;[TT22, cliques, separators, peo ] = TMFG_v2_2(r); elapsed = toc;
%     fprintf('TMFG_V2_2: %f, %f|', sum(TT22(:))+0.0, elapsed);
%     tic;[TT23, cliques, separators, peo ] = TMFG_v2_3(r); elapsed = toc;
%     fprintf('TMFG_V2_3: %f, %f\n', sum(TT23(:))+0.0, elapsed);
% end
% 
% %%%
% % RANDOM CORRELATION MATRICES TEST (RANDOM ORTHONORMAL MATRIX METHOD)
% %%%
% fprintf('Random corr matrices (Random orthogonal method)\n')
% fprintf('N|PMFG|PTIME|TMFG|TTIME|Ratio|PlanarEdges|Totaledges|Edgeratio|Avg|Stdev|Skew\n');
% for n = [50 100 150 200 300 500]
%     W = gencorr(n, 'randorth'); W = (W + W')/2; % get rid of minor differences
%     W = W .* W;
%     r = W .* W;
%     r = max(r,r');
%     fprintf('%d|', n);
%     tic;[TT, cliques, triangles, tc  ] = TMFG_mc(r); elapsed = toc;
%     fprintf('TMFG: %f, %f|', sum(TT(:))+0.0, elapsed);
%     tic;[TT1, cliques, triangles, tc  ] = TMFGT1_mc(r); elapsed = toc;
%     fprintf('TMFGT1: %f, %f|', sum(TT1(:))+0.0, elapsed);
%     tic;PP = pmfg(r);elapsed = toc;
%     fprintf('PMFG: %f, %f|', sum(PP(:))+0.0, elapsed);
%     % tic;[TT21, cliques, separators, peo ] = TMFG_v2_1(r); elapsed = toc;
%     % fprintf('TMFG_V2_1: %f, time: %f\n', sum(TT21(:))+0.0, elapsed);
%     tic;[TT22, cliques, separators, peo ] = TMFG_v2_2(r); elapsed = toc;
%     fprintf('TMFG_V2_2: %f, %f|', sum(TT22(:))+0.0, elapsed);
%     tic;[TT23, cliques, separators, peo ] = TMFG_v2_3(r); elapsed = toc;
%     fprintf('TMFG_V2_3: %f, %f\n', sum(TT23(:))+0.0, elapsed);
% end

%%%
% UNIFORM DISTRIBUTION TEST (TIMED)
%%%
% fprintf('Uniform time test\n')
% fprintf('N|TMFG|TTIME|PlanarEdges|Totaledges|Edgeratio|Avg|Stdev|Skew\n');
% for n = 100:100:10000
%     x = randraw('uniform', [0 1], [1 n*(n-1)/2]);
%     x = x ./ max(x);
%     W = triu(ones(n,n),1);
%     W(W==1) = x;
%     W = max(W, W');
%     %tic; p = pmfg(W); pte = toc;
%     % square to get TMFG weights
%     W = W .* W;
%     tic; [t c s peo] = TMFG_mc(W); tte = toc;
%     t = sqrt(t);
%     fprintf('%04d|%9.2f|%f|%d|%d|%f|', n, sum(t(:))+0, tte, 3*(n-2) , n*(n-1)/2, (3*(n-2))/(n*(n-1)/2) );
%     fprintf('%f|%f|%f\n', mean(x), std(x), skewness(x));
% end