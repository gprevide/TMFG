function [P, cliques, separators, peo] = TMFGT2_K4(W)
% Optimised algorithm, uses T2 moves only, builds a random apollonian network
% (chordal) and optimises the results by performing a permutation of the
% vertices of every 4-simplex at each stage.
% INPUT:
%   - W: a weighted network
% OUTPUT:
%   - P: the filtered TMFG graph
%   - cliques: the list of 4-cliques
%   - separators: the list of 3-cliques that are clique separators and also
%     constitute a basis for the cycle space.
%   - peo: perfect elimination ordering
% 

INITIAL_TRIANGLES = 4;

n = size(W,1); % number of vertices
P = sparse(n,n); % sparse output matrix
max_clique_gains = zeros(3*n - 6, 1); % holds the maximum gain for a triangle
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
for t = 1:INITIAL_TRIANGLES
    [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
end

for i = (INITIAL_TRIANGLES+1):n
    
    assert(size(cliques,1)== i-4);
    assert(size(separators,1) == i-5);
    
    % Get maximum local gain, nt (resp. nv) is the triangle (resp. vertex)
    % that attains the maximum gain in this configuration
    [this_gain, nt] = max(max_clique_gains);
    nv = best_vertex(nt);
    
    % Add vertex to perfect elimination order 
    peo(end + 1) = nv;
    
    % Identify the clique elements to optimize
    % This is the original clique. We want to find a permutation of the
    % clique that optimizes the sum of the weights.
    this_clique = [triangles(nt,:) nv];
    
    % Find locally optimal orientation of the new clique
    
    % First of all let's find the vertices in P that are connected to the
    % clique, excluding the vertices belonging to the triangle being
    % extended.
    tmp = find(P(triangles(nt,1), :) ~= 0);
    tmp = tmp(tmp ~= triangles(nt,1) & tmp ~= triangles(nt,2) & tmp ~= triangles(nt,3));
    neighbours_1 = tmp; % neighbours of vertex 1
    tmp = find(P(triangles(nt,2), :) ~= 0);
    tmp = tmp(tmp ~= triangles(nt,1) & tmp ~= triangles(nt,2) & tmp ~= triangles(nt,3));
    neighbours_2 = tmp; % neighbours of vertex 2
    tmp = find(P(triangles(nt,3), :) ~= 0);
    tmp = tmp(tmp ~= triangles(nt,1) & tmp ~= triangles(nt,2) & tmp ~= triangles(nt,3));
    neighbours_3 = tmp; % neighbours of vertex 3
    
    % Now let's set to zero the values in the filtered matrix, next we
    % identify the best combination and restore the links.
    P(neighbours_1, triangles(nt,1)) = 0;
    P(neighbours_2, triangles(nt,2)) = 0;
    P(neighbours_3, triangles(nt,3)) = 0;
    P(triangles(nt,1), neighbours_1) = 0;
    P(triangles(nt,2), neighbours_2) = 0;
    P(triangles(nt,3), neighbours_3) = 0;
    P(triangles(nt,:), triangles(nt,:)) = 0;
    
    local_max = 0.0;
    newsep = triangles(nt,:).';
    last_added_v = nv;
    
    for perm_clique = perms(this_clique).'
        % Assume that the last vertex is the internal one and the external
        % triangle vertices are the first three.
        int_vertex = perm_clique(4);
        ext_triangle = perm_clique(1:3);
        running_max = sum(W(neighbours_1, ext_triangle(1))) + ...
                      sum(W(neighbours_2, ext_triangle(2))) + ...
                      sum(W(neighbours_3, ext_triangle(3)));
        if running_max > local_max
            local_max = running_max;
            nv = int_vertex;
            newsep = ext_triangle;
        end
    end
    
    cliques(end + 1, :) = [newsep' nv];
    
    % Update the filtered matrix
    P(neighbours_1, newsep(1)) = W(neighbours_1, newsep(1)); 
    P(newsep(1), neighbours_1) = W(newsep(1), neighbours_1); 
    P(neighbours_2, newsep(2)) = W(neighbours_2, newsep(2)); 
    P(newsep(2), neighbours_2) = W(newsep(2), neighbours_2); 
    P(neighbours_3, newsep(3)) = W(neighbours_3, newsep(3)); 
    P(newsep(3), neighbours_3) = W(newsep(3), neighbours_3);
    P([nv; newsep], [nv; newsep]) = W([nv; newsep], [nv; newsep]); 

    if nnz(P)/2 ~= (3*i - 6)
        fprintf('Wrong update\n'); return;
    end
    test = boyer_myrvold_planarity_test(P);
    if test == 0
        fprintf('Broken planarity after update\n'); return;
    end

    % Add separator
    separators(end+1, :) = newsep;
    
    % replace triangles where some of the vertices were changed
    new_clique = [newsep' nv];
    triangles_updated = triangles;
    changed_triangles = [];
    for iv = 1:length(new_clique)
        if (new_clique(iv) ~= this_clique(iv))
            idx_t = find(triangles == this_clique(iv));
            triangles_updated(idx_t) = new_clique(iv);
            [tt ~] = ind2sub(size(triangles), idx_t);
            changed_triangles = unique([changed_triangles' tt'].');
        end
    end
    
    triangles = triangles_updated;
    
    % Replace the triangle with one of the new three, and
    triangles(nt, :) = [newsep(1) newsep(2) nv]; %#ok<*AGROW>
    % add two new triangles
    triangles(end+1, :) = [newsep(1) newsep(3) nv];
    triangles(end+1, :) = [newsep(2) newsep(3) nv];
    
    
   
    % clean cache of used values
    vertex_list = setdiff(vertex_list, last_added_v);
    
    %
    % Update max_gains reflecting the fact that nv cannot be used any
    % longer. Find other triangles coupled with nv. vertex_list does not
    % contain nv at this stage.
    % 
    if ~isempty(vertex_list)
        for t = find(best_vertex == last_added_v).'
            [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
        end
    end
    
    %
    % Now update the gains list with the new triangles 
    %
    max_clique_gains(nt) = 0.0;
    best_vertex(nt) = 0;
    ct = size(triangles, 1);
    if ~isempty(vertex_list)
        for t = [nt (ct-1) ct]
            [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
        end
    end
    
    % and the changed triangles
    if ~isempty(vertex_list)
        for t = changed_triangles'
            if max_clique_gains(t) ~= 0
                [max_clique_gains(t) best_vertex(t)] = get_best_gain(vertex_list, triangles(t,:), W);
            end
        end
    end
end
% %---------------
% %
% % Apply a last round of optimisation by flipping triangles
% %
% tot_tr = size(triangles, 1);
% prev_sum = sum(P(:)) + 0.0;
% go = 1;
% while go
%     go = 0;
%     for cur_t = 1:tot_tr
%         for oth_t = (cur_t+1):tot_tr
%             if is_neighbour(triangles, cur_t, oth_t)
%                 [triangles P flipped] = flip(cur_t, oth_t, triangles, W, P);
%             end
%         end
%     end
%     new_sum = sum(P(:)) + 0.0;
%     if new_sum > prev_sum + 0.00000001
%         go = 1; 
%         %fprintf('\nT1 >>> Previous %10.5f Now %10.5f nnz %d planar %d\n', (prev_sum+0.0)/2, (new_sum + 0.0) /2,(nnz(P)+0.0)/2,boyer_myrvold_planarity_test(sparse(P)));
%         prev_sum = new_sum;
%     end
% end
% %---------------
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

function [triangles_out P_out flipped] = flip(t1, t2, triangles, W, P)
    flipped = false;
    C = intersect(triangles(t1, :), triangles(t2, :));
    assert(length(C)==2);
    v1 = C(1); v3 = C(2);
    v2 = setdiff(triangles(t1, :), [v1 v3]);
    v4 = setdiff(triangles(t2, :), [v1 v3]);
    if (P(v2, v4) + 0) == 0 && (W(v2, v4) > W(v1, v3))
        triangles_out = triangles;
        triangles_out(t1, :) = [v1 v2 v4];
        triangles_out(t2, :) = [v2 v3 v4];
        P_out = P;
        P_out(v2, v4) = W(v2, v4);
        P_out(v4, v2) = W(v4, v2);
        P_out(v1, v3) = 0;
        P_out(v3, v1) = 0;
        flipped = true;
        %fprintf('flipped %d %d\n', t1, t2);
    else
        P_out = P;
        triangles_out = triangles;
    end
end

function isn = is_neighbour(triangles, t1, t2)
    cnt = length(intersect(triangles(t1, :), triangles(t2, :)));
    if cnt == 2 ; isn = true; else isn = 0 ;end
end

