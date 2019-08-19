function [ P, cliques, triangles, peo ] = TMFG_mc(W)
% BASE Algorithm, uses only T2, builds a random apollonian network
% (chordal)
% INPUT:
%   - W: a weighted network
% OUTPUT:
%   - P: the filtered TMFG graph
%   - cliques: the list of 4-cliques
%   - separators: the list of 3-cliques that are clique separators and also
%     constitute a basis for the cycle space.
%   - junction_tree: one possible junction tree for the graph
%   - peo: perfect elimination ordering

n = size(W,1); % number of vertices
P = sparse(n,n); % sparse matrix
max_clique_gains = zeros(3*n - 6, 1);
best_vertex = zeros(3*n - 6, 1);

cliques = [];
separators = [];
triangles = [];

% Get first simplex and populate the first 4 triangles in the basis
cliques(1, :) = max_clique(W);
vertex_list = setdiff(1:n, cliques(1,:));
triangles(1,:) = cliques(1, [1 2 3]); 
triangles(2,:) = cliques(1, [1 2 4]); 
triangles(3,:) = cliques(1, [1 3 4]); 
triangles(4,:) = cliques(1, [2 3 4]); 

peo = cliques(1, :); % start perfect elimination order with 1st clique

W(1:(n+1):n^2) = 0;

P(peo, peo) = W(peo, peo);

% init gain matrix
for t = 1:4
    [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
end

for i = 1:(n-4)
    % get maximum local gain
    [~, nt] = max(max_clique_gains);
    nv = best_vertex(nt);
    %disp(nv);
    
    peo(end + 1) = nv;
    % add clique
    cliques(end+1, :) = [nv triangles(nt,:)];
    % add separators
    newsep = triangles(nt, :);
    P([nv newsep], [nv newsep]) = W([nv newsep], [nv newsep]);
%     if ~boyer_myrvold_planarity_test(sparse(P)) 
%         fprintf('n: %d\n', i);
%         break;
%     end
    separators(end+1, :) = newsep;
    % replace triangles
    triangles(nt, :) = [newsep(1) newsep(2) nv];
    % add two new triangles
    triangles(end+1, :) = [newsep(1) newsep(3) nv];
    triangles(end+1, :) = [newsep(2) newsep(3) nv];
    % clean cache of used values
    vertex_list = setdiff(vertex_list, nv);
    % update max gains where the vertex nv was involved
    % 
    if length(vertex_list) > 0
        for t = find(best_vertex == nv).'
            [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
        end
    end
    max_clique_gains(nt) = 0;
    ct = size(triangles, 1);
    if length(vertex_list) > 0
        for t = [nt (ct-1) ct]
            [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
        end
    end
end

end

function [gain vertex] = get_best_gain(vertex_list, triangle, W)
    gvec(vertex_list) = W(vertex_list, triangle(1)) + W(vertex_list, triangle(2)) + W(vertex_list, triangle(3));
    [gain vertex] = max(gvec);
end

function gvec = gain_vec(vertex_list, triangle, W)
    gvec = W(vertex_list, triangle(1)) + W(vertex_list, triangle(2)) + W(vertex_list, triangle(3));
end

function cl = max_clique(W)
    v = sum(W.*(W>mean(W(:))),2);
    % v = sum(W .* W);
    [~, sortindex] = sort(v, 'descend');
    cl = sortindex(1:4);
end