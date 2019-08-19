function [ P, cliques, triangles, tc  ] = TMFGT1_mc(W)

n = size(W,1);   % number of vertices
P = sparse(n,n); % sparse matrix
max_clique_gains = zeros(3*n - 6, 1);
best_vertex      = zeros(3*n - 6, 1);
tc = sparse(n,n); % holds the contact structure between triangles (only in basis)

cliques = [];
triangles = [];

% Get first simplex and populate the first 4 triangles in the basis
cliques(1, :) = max_clique(W);
vertex_list = setdiff(1:n, cliques(1,:));
triangles(1,:) = cliques(1, [1 2 3]); 
triangles(2,:) = cliques(1, [1 2 4]); 
triangles(3,:) = cliques(1, [1 3 4]); 
triangles(4,:) = cliques(1, [2 3 4]); 

% Populate the contact structure 
tc([1 2 3 4], [1 2 3 4]) = 1;
tc(1:(n+1):n^2) = 0;

peo = cliques(1, :); 

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
    
   
    % add clique
    cliques(end+1, :) = [nv triangles(nt,:)];
    % add separators
    newsep = triangles(nt, :);
    P([nv newsep], [nv newsep]) = W([nv newsep], [nv newsep]);
%     if ~boyer_myrvold_planarity_test(sparse(P)) 
%         fprintf('n: %d\n', i);
%         break;
%     end
    % replace triangle
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
    possible_neighbours = find(tc(nt,:)); % find neighbours of triangle going out of the basis
    tc(nt, :) = 0; tc(:,nt) = 0; % remove contact structure for old triangle
    ct = size(triangles, 1);
    if length(vertex_list) > 0
        for t = [nt (ct-1) ct]
            [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
        end
    end
    new_triangles = [nt (ct-1) ct];
    tc(new_triangles, new_triangles) = 1;
    tc(new_triangles(1), new_triangles(1)) = 0;
    tc(new_triangles(2), new_triangles(2)) = 0;
    tc(new_triangles(3), new_triangles(3)) = 0;
    for t1 = possible_neighbours
        for t2 = new_triangles
            tc(t1, t2) = is_neighbour(triangles, t1, t2);
            tc(t2, t1) = tc(t1, t2);
            if tc(t1, t2)
                [triangles P flipped] =  flip(t1, t2, triangles, W, P);
            end
        end
    end
   
    if ~isempty(vertex_list)
        for t = new_triangles
            [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
        end
    end

end

end

function [gain vertex] = get_best_gain(vertex_list, triangle, W)
    gvec(vertex_list) = W(vertex_list, triangle(1)) + W(vertex_list, triangle(2)) + W(vertex_list, triangle(3));
    [gain vertex] = max(gvec);
end

function cl = max_clique(W)
    v = sum(W.*(W>mean(W(:))),2);
    % v = sum(W .* W);
    [~, sortindex] = sort(v, 'descend');
    cl = sortindex(1:4);
end

function [triangles_out P_out flipped] = flip(t1, t2, triangles, W, P)
    flipped = false;
    C = intersect(triangles(t1, :), triangles(t2, :));
    v1 = C(1); v3 = C(2);
    v2 = setdiff(triangles(t1, :), [v1 v3]);
    v4 = setdiff(triangles(t2, :), [v1 v3]);
%     fprintf('%f %d %d\n',P(v2,v4)+0, v2, v4);
    if (~isempty(v2) && ~isempty(v4) && ...
        ~isempty(v1) && ~isempty(v3) && ...
        (P(v2, v4) + 0) == 0 && (W(v2, v4) > W(v1, v3)))
        triangles_out = triangles;
        triangles_out(t1, :) = [v1 v2 v4];
        triangles_out(t2, :) = [v2 v3 v4];
        P_out = P;
        P_out(v2, v4) = W(v2, v4);
        P_out(v4, v2) = W(v4, v2);
        P_out(v1, v3) = 0;
        P_out(v3, v1) = 0;
        flipped = true;
    else
        P_out = P;
        triangles_out = triangles;
    end
end

function isn = is_neighbour(triangles, t1, t2)
    cnt = length(intersect(triangles(t1, :), triangles(t2, :)));
    if cnt == 2 ; isn = true; else isn = 0 ;end
end

