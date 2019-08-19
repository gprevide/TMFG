addpath('C:\Users\Guido\Documents\Bauhaus\PhD\matlab_bgl', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\pmfg', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\BCT', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\RANDRAW')
    
fileid = fopen('.\Paper_test_Uniform.txt', 'a');
%----------------------|
fprintf(fileid, 'Uniform [0,1]\n');
fprintf(fileid, 'iter|TMFG_Val|TMFG_elapsed|TMFGT1_Val|TMFGT1_elapsed|TMFGT2_K4_Val|TMFGT2_K4_elapsed|PMFG_Val|PMFG_elapsed\n');
for iter = 1:100 % 100 simulations
    for n = [400] % of 400x400 matrices
        x = randraw('uniform', [0 1], [1 n*(n-1)/2]);
        x = (x - min(x)) ./ (max(x) - min(x));
        W = triu(ones(n,n),1);
        W(W==1) = x;
        r = W .* W;
        r = max(r,r');
        fprintf(fileid, '%d|', iter);
        tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
        fprintf(fileid, '%f|%f|', sum(TT(:))+0.0, elapsed);
        tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
        fprintf(fileid,'%f|%f|', sum(TT1(:))+0.0, elapsed);
        tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
        fprintf(fileid, '%f|%f|', sum(TT22(:))+0.0, elapsed);
        tic;PP = pmfg(r);elapsed = toc;
        fprintf(fileid,'%f|%f\n', sum(PP(:))+0.0, elapsed);
    end
end
fclose(fileid);

%%%
% Beta
%%%

fileid = fopen('.\Paper_test_Beta.txt', 'a');
for params = [0.5 10; 10 0.5]
    %----------------------|
    fprintf(fileid, 'Beta alpha %f beta %f\n', params(1), params(2));
    fprintf(fileid, 'iter|TMFG_Val|TMFG_elapsed|TMFGT1_Val|TMFGT1_elapsed|TMFGT2_K4_Val|TMFGT2_K4_elapsed|PMFG_Val|PMFG_elapsed\n');
    for iter = 1:100 % 100 simulations
        for n = [400] % of 400x400 matrices
            x = randraw('beta', [params(1) params(2)], [1 n*(n-1)/2]);
            x = (x - min(x)) ./ (max(x) - min(x));
            W = triu(ones(n,n),1);
            W(W==1) = x;
            r = W .* W;
            r = max(r,r');
            fprintf(fileid, '%d|', iter);
            tic;[TT, cliques, triangles, tc  ] = TMFG(r); elapsed = toc;
            fprintf(fileid, '%f|%f|', sum(TT(:))+0.0, elapsed);
            tic;[TT1, cliques, triangles, tc  ] = TMFGT1(r); elapsed = toc;
            fprintf(fileid,'%f|%f|', sum(TT1(:))+0.0, elapsed);
            tic;[TT22, cliques, separators, peo ] = TMFGT2_K4(r); elapsed = toc;
            fprintf(fileid, '%f|%f|', sum(TT22(:))+0.0, elapsed);
            tic;PP = pmfg(r);elapsed = toc;
            fprintf(fileid,'%f|%f\n', sum(PP(:))+0.0, elapsed);
        end
    end
end
fclose(fileid);

