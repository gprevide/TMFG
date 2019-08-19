function [ corrmat ] = gencorr( n, varargin )
%GENCORR Summary of this function goes here
%   Detailed explanation goes here
corrmat = [];

if nargin > 1
    method = varargin(1);
else
    disp('User must specify a method');
    return;
end

if strcmp(method, 'factors')   
    if nargin > 2
        sys_fac_no = cell2mat(varargin(2));
    else
        sys_fac_no = ceil(n/10);
    end
    if nargin > 3
        df = cell2mat(varargin(3));
    else
        df = ceil(n * 10);
    end
    X = zeros(n, df);
    sys_factors = randn(sys_fac_no, df);
    for a = 1:n
        rsq = rand(); % r-squared
        weights = -1 + 2*rand(sys_fac_no,1);
        weights =  sqrt(rsq) * weights / sqrt(sum (weights .* weights));
        X(a, :) = weights' * sys_factors  + sqrt(1-rsq)*rand(1, df);
    end
    corrmat = corrcoef(X');
end

if strcmp(method, 'randorth')
    Q = orth(rand(n));
    eigenvalues = rand(n,1); eigenvalues * n / sum(eigenvalues);
    corrmat = Q * diag(eigenvalues) * Q';
end

end

