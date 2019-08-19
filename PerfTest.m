% 28 Mar 2015
% Performance and timing test for TMFG, PMFG, etc...
% ---

addpath('C:\Users\Guido\Documents\Bauhaus\PhD\matlab_bgl', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\pmfg', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\BCT', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\RANDRAW');
    
%fileid = 1 ;%fopen('.\Paper_test_all.txt', 'a');
fileid = fopen('.\Paper_test_New.txt', 'a');
samples = 1;
sizes = [200];

rng(1234);

%%%
% POWER DISTRIBUTION TEST (PARETO)
%%%
for alpha = [1 2]
    %----------------------|
    fprintf(fileid, 'Pareto Alpha %d\n', alpha);
    fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
    for iter = 1:samples % e.g. 100 simulations
        for n = sizes % e.g. 400x400 matrices
            x = randraw('pareto', [1, alpha], [1 n*(n-1)/2]);
            x = (x - min(x)) ./ (max(x) - min(x));
            W = triu(ones(n,n),1);
            W(W==1) = x;
            r = W .* W;
            r = max(r,r');
            fprintf(fileid, '%d,', n);
            filename = sprintf('pareto_%d_%d.csv', alpha, n);
            dlmwrite(filename,r);
            tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
            tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
            tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
            tic;PP = pmfg(r);elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
            r = tril(r, -1); v = sort(r(:), 'descend'); upperbound = 2 * sum(v(1:3*(n-2)));
            fprintf(fileid,'%f\n', upperbound);
        end
    end
end

%%%
% BETA DISTRIBUTION TEST
%%%
for  params = [0.5 3; 3 0.5]
    %----------------------|
    fprintf(fileid, 'Beta alpha %f beta %f\n', params(1), params(2));
    fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
    for iter = 1:samples % e.g. 100 simulations
        for n = sizes % e.g. 400x400 matrices
            x = randraw('beta', [params(1) params(2)], [1 n*(n-1)/2]);
            x = (x - min(x)) ./ (max(x) - min(x));
            W = triu(ones(n,n),1);
            W(W==1) = x;
            r = W .* W;
            r = max(r,r');
            fprintf(fileid, '%d,', n);
            filename = sprintf('beta_%f_%f_%d.csv', params(1), params(2), n);
            dlmwrite(filename,r);
            tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
            tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
            tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
            tic;PP = pmfg(r);elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
            r = tril(r, -1); v = sort(r(:), 'descend'); upperbound = 2 * sum(v(1:3*(n-2)));
            fprintf(fileid,'%f\n', upperbound);
        end
    end
end

%%%
% UNIFORM DISTRIBUTION TEST
%%%
fprintf(fileid, 'Uniform [0,1]\n');
fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
for iter = 1:samples % e.g. 100 simulations
    for n = sizes % e.g. 400x400 matrices
        x = randraw('uniform', [0 1], [1 n*(n-1)/2]);
        x = (x - min(x)) ./ (max(x) - min(x));
        W = triu(ones(n,n),1);
        W(W==1) = x;
        r = W .* W;
        r = max(r,r');
        fprintf(fileid, '%d,', n);
        filename = sprintf('uniform_%d.csv', n);
        dlmwrite(filename,r);
        tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
        fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
        tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
        fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
        tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
        fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
        tic;PP = pmfg(r);elapsed = toc;
        fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
        r = tril(r, -1); v = sort(r(:), 'descend'); upperbound = 2 * sum(v(1:3*(n-2)));
        fprintf(fileid,'%f\n', upperbound);
    end
end

%%%
% FACTOR DISTRIBUTION TEST
%%%
for factors = [20 50 100]
    fprintf(fileid, 'Corr with %d factors\n', factors);
    fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
    for iter = 1:samples % e.g. 100 simulations
        for n = sizes % e.g. 400x400 matrices
            W = gencorr(n, 'factors', factors); W = (W + W')/2; % get rid of minor differences
            r = W .* W;
            r = max(r,r');
            fprintf(fileid, '%d,', n);
            filename = sprintf('factors_%d_%d.csv', factors, n);
            dlmwrite(filename,r);
            tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
            tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
            tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
            tic;PP = pmfg(r);elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
            r = tril(r, -1); v = sort(r(:), 'descend'); upperbound = 2 * sum(v(1:3*(n-2)));
            fprintf(fileid,'%f\n', upperbound);
        end
    end
end


%%%
% REAL MATRIX
%%%
% load('342USstocks.mat')
% % there is one crazy price! YRC WORLDWIDE INC i=296 (but seems correct, it is an effect of the adjusting...)
% Pr = Prices(:,[1:295,297:end]);
% wi = 1000;
% r = corrcoef(diff(log(Pr(1:wi,:))));
% k = 0;
% 
% t=wi;
% k=1;
% % from here can restart any time
% t0=t;
% k=k-1;
% siz = size(Prices,2);
% 
% fprintf(fileid, 'Real correlation matrix sample\n');
% fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
% for t=t0:20:size(Pr,1)
%     k = k+1;
%     t_d(k) = dates(t);
%     r = corrcoef(diff(log(Pr((t-wi+1):t, :))));
%     R = r.^2; R=(R+R')/2;
%     fprintf(fileid, '%d,', siz);
%     tic;[TT, cliquesT, trianglesT, tcT  ] = TMFG(R); elapsed = toc; % re-buid TMFG T2 all times
%     fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
%     tic;[TT1, cliquesT1, trianglesT1, tcT1  ] = TMFGT1(R); elapsed = toc;% re-buid TMFG T1 all times
%     fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
%     tic;[TT22, ~, ~, ~ ] = TMFGT2_K4(R); elapsed = toc; %last version of code where swaps are done while building the graph
%     fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
%     tic;PP = pmfg(R); elapsed = toc;% re-buid PMFG all time
%     fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
%     R = tril(R, -1); v = sort(R(:), 'descend'); upperbound = 2 * sum(v(1:3*(siz-2)));
%     fprintf(fileid,'%f\n', upperbound);
% end      
% 
% %%%
% % POWER DISTRIBUTION TEST (PARETO) INCREASING SIZE
% %%%
% for alpha = [1]
%     %----------------------|
%     fprintf(fileid, 'Pareto Alpha %d\n', alpha);
%     fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
%     for iter = 1:samples % e.g. 100 simulations
%         for n = [50 100 150 300 500 700] % e.g. 400x400 matrices
%             x = randraw('pareto', [1, alpha], [1 n*(n-1)/2]);
%             x = (x - min(x)) ./ (max(x) - min(x));
%             W = triu(ones(n,n),1);
%             W(W==1) = x;
%             r = W .* W;
%             r = max(r,r');
%             fprintf(fileid, '%d,', n);
%             tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
%             fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
%             tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
%             fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
%             tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
%             fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
%             tic;PP = pmfg(r);elapsed = toc;
%             fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
%             r = tril(r, -1); v = sort(r(:), 'descend'); upperbound = 2 * sum(v(1:3*(n-2)));
%             fprintf(fileid,'%f\n', upperbound);
%         end
%     end
% end

%%%
% BETA DISTRIBUTION TEST INCREASING SIZE
%%%
for  params = [0.5 ; 3 ]
    %----------------------|
    fprintf(fileid, 'Beta alpha %f beta %f\n', params(1), params(2));
    fprintf(fileid, 'size,TMFG_Val,TMFG_elapsed,TMFGT1_Val,TMFGT1_elapsed,TMFGT2_K4_Val,TMFGT2_K4_elapsed,PMFG_Val,PMFG_elapsed,UpperBound\n');
    for iter = 1:samples % e.g. 100 simulations
        for n = n = [50 100 150 300 500 700]  % e.g. 400x400 matrices
            x = randraw('beta', [params(1) params(2)], [1 n*(n-1)/2]);
            x = (x - min(x)) ./ (max(x) - min(x));
            W = triu(ones(n,n),1);
            W(W==1) = x;
            r = W .* W;
            r = max(r,r');
            fprintf(fileid, '%d,', n);
            tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT(:))+0.0, elapsed);
            tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(TT1(:))+0.0, elapsed);
            tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
            fprintf(fileid, '%f,%f,', sum(TT22(:))+0.0, elapsed);
            tic;PP = pmfg(r);elapsed = toc;
            fprintf(fileid,'%f,%f,', sum(PP(:))+0.0, elapsed);
            r = tril(r, -1); v = sort(r(:), 'descend'); upperbound = 2 * sum(v(1:3*(n-2)));
            fprintf(fileid,'%f\n', upperbound);
        end
    end
end

if fileid ~=1, fclose(fileid), end;