function [ P, triangles, tc  ] = TMFGLoz(w, num_swaps)

n = size(w,1);   % number of vertices
P = sparse(n,n); % sparse matrix

max_loz_gains = sparse(3*n,3*n);
best_vertex   = sparse(3*n,3*n);
best_op       = sparse(3*n,3*n);

tc = sparse(false(3*n,3*n)); % holds the contact structure between triangles (only in basis)
triangles =  sparse(false(3*n, n));
triangles_changed = [];

% Get first simplex and populate the first 4 triangles in the basis
K_4 = max_clique(w);
vertex_list = setdiff(1:n, K_4);
triangles(K_4([1 2 3]), 1) = true; 
triangles(K_4([1 2 4]), 2) = true; 
triangles(K_4([1 3 4]), 3) = true; 
triangles(K_4([2 3 4]), 4) = true; 
triangles_changed([1 2 3 4]) = 1;

% Populate the contact structure 
tc([1 2 3 4], [1 2 3 4]) = 1;
tc = tril(tc, -1);

w(1:(n+1):n^2) = 0;

P(K_4, K_4) = w(K_4, K_4);

% init gain matrix
% Lower triangular
% Loop over adjacent triangles and find best operation and best gain
[row col val] = find(tril(tc));

for idx = 1:numel(row)
    [max_loz_gains(row(idx), col(idx)) best_vertex(row(idx), col(idx)) best_op(row(idx), col(idx))] = ...
            get_loz_gains(w, triangles, row(idx), col(idx), vertex_list, P);
end

for i = 1:(n-4)
    % get maximum local gain
    [~, nt] = max(max_loz_gains(:)); % this is a matrix maximum
    nv = best_vertex(nt); % nv is the vertex to insert
    op = best_op(nt); % op is the operation to apply
    [t1 t2] = ind2sub(size(max_loz_gains), nt);
    %fprintf('op %d v %d\n', op + 0, nv + 0);
    [P triangles tc triangles_changed] = apply_op(op, w, nv, P, t1, t2, triangles, tc);
    vertex_list = setdiff(vertex_list, nv);
    % delete gains where deleted lozenges were involved
    max_loz_gains = max_loz_gains .* (tc ~= 0);
    best_vertex = best_vertex .* (tc ~= 0);
    best_op = best_op .* (tc ~= 0);
    
    % update max gains where the vertex nv was involved
    if ~isempty(vertex_list)
        [row col val] = find(best_vertex == nv);
        for idx = 1:numel(row)
            max_loz_gains(row(idx), col(idx))= 0;
            best_vertex(row(idx), col(idx))= 0;
            best_op(row(idx), col(idx))= 0;
            if tc(row(idx), col(idx)) == true
                [max_loz_gains(row(idx), col(idx)) best_vertex(row(idx), col(idx)) best_op(row(idx), col(idx))] = ...
                    get_loz_gains(w, triangles, row(idx), col(idx), vertex_list, P);
            end
        end
        for it = 1:num_swaps
            [P triangles tc triangles_changed] = apply_swaps(w, P, triangles, tc);
        end
        [row col val] = find(tril(tc));
        for idx = 1:numel(row)
            if triangles_changed(row(idx)) || triangles_changed(col(idx)) || max_loz_gains(row(idx), col(idx)) == 0
                [max_loz_gains(row(idx), col(idx)) best_vertex(row(idx), col(idx)) best_op(row(idx), col(idx))] = ...
                    get_loz_gains(w, triangles, row(idx), col(idx), vertex_list, P);
                triangles_changed(row(idx)) = 0;
                triangles_changed(col(idx)) = 0;
            end
        end
    end
end

end

function cl = max_clique(W)
    v = sum(W.*(W>mean(W(:))),2);
    % v = sum(W .* W);
    [~, sortindex] = sort(v, 'descend');
    cl = sortindex(1:4);
end

function isn = is_neighbour(triangles, t1, t2)
    cnt = length(find(triangles(:, t1) & triangles(:, t2)));
    if cnt == 2 ; isn = true; else isn = 0 ;end
end

function [loz_gain vertex op] =  get_loz_gains(w, triangles, t1, t2, vertex_list, P)
    gvec = zeros(numel(vertex_list), 6);
    % ----
    N = tsetdiff(triangles, t1, t2);
    S = tsetdiff(triangles, t2, t1);
%     NS = txor(triangles, t1, t2); 
%     N = NS(1);
%     S = NS(2);
    tmp = tintersect(triangles, t1, t2);
    if (numel(tmp) ~= 2)
        fprintf('numel %d\n', numel(tmp));
    end
    W = tmp(1);
    E = tmp(2);
    % ----
    % operation 1
    gvec(vertex_list, 1) = w(vertex_list, N) + w(vertex_list, W) + w(vertex_list, E);
    % operation 2
    gvec(vertex_list, 2) = w(vertex_list, N) + w(vertex_list, S) + w(vertex_list, E) + w(vertex_list, W) - w(W,E);
    % operation 3
    gvec(vertex_list, 3) = w(vertex_list, S) + w(vertex_list, W) + w(vertex_list, E);
    if (P(N,S) == 0)
        % operation 4
        gvec(vertex_list, 4) = w(vertex_list, N) + w(vertex_list, S) + w(vertex_list, W) + w(N,S) - w(W,E);
        % operation 5
        gvec(vertex_list, 5) = w(vertex_list, N) + w(vertex_list, S) + w(vertex_list, E) + w(N,S) - w(W,E);
    else
        % operation 4
        gvec(vertex_list, 4) = 0;
        % operation 5
        gvec(vertex_list, 5) = 0;
    end
    
    idx = find(gvec == max(max(gvec)));
    loz_gain = gvec(idx);
    [vertex op] = ind2sub(size(gvec), idx);
end

function [P triangles tc triangles_changed] = apply_op(op, w, nv, P, t1, t2, triangles, tc)
    triangles_changed = [];
    N = tsetdiff(triangles, t1, t2);
    S = tsetdiff(triangles, t2, t1);
    tmp = tintersect(triangles, t1, t2);
    if (numel(tmp) <2)
        fprintf('hhh\n')
    end
    W = tmp(1);
    E = tmp(2);
    C = nv;
    switch op
        case 1
            % adjust triangles
            triangles( :, t1) = false;
            triangles([W E C], t1) = true;
            triangles_changed(t1) = 1;
            [~, y] = find(triangles) ; t3 = max(y) +1;
            t4 = t3 + 1;
            triangles( :, t3) = false;
            triangles([W C N], t3) = true;
            triangles_changed(t3) = 1;
            triangles( :, t4) = false;
            triangles([E C N], t4) = true;
            triangles_changed(t4) = 1;
            % adjust tc externally
            %   find triangles bordering with t1
            for tn = find_neib(tc, t1)
                if tn == t2
                    % do nothing because t1 and t2 are neighbours anyway by
                    % construction
                end
                if is_neighbour(triangles, tn, t3)
                    tc = connecttril(tc, tn, t3);
                    tc(max(tn, t1), min(tn,t1)) = 0;
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t4)
                    tc = connecttril(tc, tn, t4);
                    tc(max(tn, t1), min(tn,t1)) = 0;
                    triangles_changed(tn) = 1;
                end  
            end
            % adjust tc internally
            tc(t1,:) =0; tc(:,t1)= 0;
            tc = connecttril(tc, t1, t2);
            tc = connecttril(tc, t1, t3);
            tc = connecttril(tc, t1, t4);
            tc = connecttril(tc, t3, t4);
            % adjust P
            P(W, C) = w(W, C); P(C, W) = w(C, W);
            P(E, C) = w(E, C); P(C, E) = w(C, E);
            P(N, C) = w(N, C); P(C, N) = w(C, N);
        case 2
            triangles( :, t1) = false;
            triangles([W C N], t1) = true;
            triangles_changed(t1) = 1;
            triangles( :, t2) = false;
            triangles([W C S], t2) = true;
            triangles_changed(t2) = 1;
            [~, y] = find(triangles) ; t3 = max(y) +1;
            t4 = t3 + 1;
            triangles( :, t3) = false;
            triangles([N C E], t3) = true;
            triangles_changed(t3) = 1;
            triangles( :, t4) = false;
            triangles([S C E], t4) = true;
            triangles_changed(t4) = 1;
            % adjust tc externally
            % boundary triangles
            tb1 =  find_neib(tc, t1);
            tb2 =  find_neib(tc, t2);
            tb = union(tb1, tb2);
            tc(t1, :) = 0; tc(:, t1) = 0;
            tc(t2, :) = 0; tc(:, t2) = 0;
            for tn = tb
                if is_neighbour(triangles, tn, t1)
                    tc(max(tn, t1), min(tn,t1)) = 1;
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t2)
                    tc(max(tn, t2), min(tn,t2)) = 1;
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t3)
                    tc(max(tn, t3), min(tn,t3)) = 1;
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t4)
                    tc(max(tn, t4), min(tn,t4)) = 1;
                    triangles_changed(tn) = 1;
                end
            end
            % adjust tc internally
            tc = connecttril(tc, t1, t3);
            tc = connecttril(tc, t1, t2);
            tc = connecttril(tc, t2, t4);
            tc = connecttril(tc, t3, t4);
            % adjust P
            P(W, E) = 0; P(E, W) = 0;
            P(W, C) = w(W, C); P(C, W) = w(C, W);
            P(E, C) = w(E, C); P(C, E) = w(C, E);
            P(S, C) = w(S, C); P(C, S) = w(C, S);
            P(N, C) = w(N, C); P(C, N) = w(C, N);
        case 3
            % adjust triangles
            triangles( :, t2) = false;
            triangles([W E C], t2) = true;
            triangles_changed(t2) = 1;
            [~, y] = find(triangles) ; t3 = max(y) +1;
            t4 = t3 + 1;
            triangles( :, t3) = false;
            triangles([W C S], t3) = true;
            triangles_changed(t3) = 1;
            triangles( :, t4) = false;
            triangles([E C S], t4) = true;
            triangles_changed(t4) = 1;
            % adjust tc externally
            %   find triangles bordering with t1
            for tn = find_neib(tc, t2)
                if tn == t1
                    % do nothing
                end
                if is_neighbour(triangles, tn, t3)
                    tc(max(tn, t3), min(tn,t3)) = 1;
                    tc(max(tn, t2), min(tn,t2)) = 0;
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t4)
                    tc(max(tn, t4), min(tn,t4)) = 1;
                    tc(max(tn, t2), min(tn,t2)) = 0;
                    triangles_changed(tn) = 1;
                end  
            end
            % adjust tc internally
            tc(t2,:) =0; tc(:,t2)= 0;
            tc = connecttril(tc, t1, t2);
            tc = connecttril(tc, t2, t3);
            tc = connecttril(tc, t2, t4);
            tc = connecttril(tc, t3, t4);
            % adjust P
            P(W, C) = w(W, C);
            P(C, W) = w(C, W);
            P(E, C) = w(E, C);
            P(C, E) = w(C, E);
            P(S, C) = w(S, C);
            P(C, S) = w(C, S);
        case 4
            % t1
            triangles( :, t1) = false;
            triangles([W C N], t1) = true;
            triangles_changed(t1) = 1;
            % t2
            triangles( :, t2) = false;
            triangles([W C S], t2) = true;
            triangles_changed(t2) = 1;
            % increment triangles count
            [~, y] = find(triangles) ; t3 = max(y) +1;
            t4 = t3 + 1;
            % t3
            triangles( :, t3) = false;
            triangles([N S C], t3) = true;
            triangles_changed(t3) = 1;
            % t4
            triangles( :, t4) = false;
            triangles([N S E], t4) = true;
            triangles_changed(t4) = 1;
            % process neighboring triangles
            tb1 = find_neib(tc, t1);
            tb2 = find_neib(tc, t2);
            tb = union(tb1, tb2);
            % disconnect neighboring triangles
            tc(t1, :) = 0; tc(:, t1) = 0;
            tc(t2, :) = 0; tc(:, t2) = 0;
            % reconnect external triangles
            for tn = tb
                if is_neighbour(triangles, tn, t1)
                    tc = connecttril(tc, tn, t1);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t2)
                    tc = connecttril(tc, tn, t2);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t3)
                    tc = connecttril(tc, tn, t3);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t4)
                    tc = connecttril(tc, tn, t4);
                    triangles_changed(tn) = 1;
                end  
            end
            % adjust tc internally
            tc = connecttril(tc, t1, t2);
            tc = connecttril(tc, t1, t3);
            tc = connecttril(tc, t2, t3);
            tc = connecttril(tc, t3, t4);
            % adjust P
            P(W, E) = 0; P(E, W) = 0;
            P(N, S) = w(N, S); P(S, N) = w(S, N);
            P(S, C) = w(S, C); P(C, S) = w(C, S);
            P(N, C) = w(N, C); P(C, N) = w(C, N);
            P(W, C) = w(W, C); P(C, W) = w(C, W);
        case 5
            % t1
            triangles( :, t1) = false;
            triangles([N E C], t1) = true;
            triangles_changed(t1) = 1;
            % t2
            triangles( :, t2) = false;
            triangles([E C S], t2) = true;
            triangles_changed(t2) = 1;
            % increment triangles count
            [~, y] = find(triangles) ; t3 = max(y) +1;
            t4 = t3 + 1;
            % t3
            triangles( :, t3) = false;
            triangles([N S C], t3) = true;
            triangles_changed(t3) = 1;
            % t4
            triangles( :, t4) = false;
            triangles([N S W], t4) = true;
            triangles_changed(t4) = 1;
            % process neighboring triangles
            tb1 = find_neib(tc, t1);
            tb2 = find_neib(tc, t2);
            tb = union(tb1, tb2);
            % disconnect neighboring triangles
            tc(t1, :) = 0; tc(:, t1) = 0;
            tc(t2, :) = 0; tc(:, t2) = 0;
            % reconnect external triangles
            for tn = tb
                if is_neighbour(triangles, tn, t1)
                    tc = connecttril(tc, tn, t1);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t2)
                    tc = connecttril(tc, tn, t2);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t3)
                    tc = connecttril(tc, tn, t3);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t4)
                    tc = connecttril(tc, tn, t4);
                    triangles_changed(tn) = 1;
                end  
            end
            % adjust tc internally
            tc = connecttril(tc, t1, t2);
            tc = connecttril(tc, t1, t3);
            tc = connecttril(tc, t2, t3);
            tc = connecttril(tc, t3, t4);
            % adjust P
            P(W, E) = 0; P(E, W) = 0;
            P(N, S) = w(N, S); P(S, N) = w(S, N);
            P(S, C) = w(S, C); P(C, S) = w(C, S);
            P(N, C) = w(N, C); P(C, N) = w(C, N);
            P(E, C) = w(E, C); P(C, E) = w(C, E);
        otherwise
            fprintf('some serious error here');
            return;
    end
end

function [P triangles tc triangles_changed] = apply_swaps(w, P, triangles, tc)
    [row col] = find(tril(tc));
    for idx = 1:numel(row)
        t1 = row(idx); t2 = col(idx);
        if ~tc(t1, t2) % the contact might have changed in the loop
            continue;
        end
        N = tsetdiff(triangles, t1, t2);
        S = tsetdiff(triangles, t2, t1);
        tmp = tintersect(triangles, t1, t2);
        W = tmp(1);
        E = tmp(2);
        triangles_changed(t1) = 0; triangles_changed(t2) = 0;
        if w(N, S) > w(W, E) && P(N, S) == 0
            P(W, E) = 0; P(E, W) = 0;
            P(N, S) = w(N, S); P(S, N) = w(S, N);
            triangles(:, t1) = false;
            triangles([N W S], t1) = true;
            triangles_changed(t1) = true;
            triangles(:, t2) = false;
            triangles([N S E], t2) = true;
            triangles_changed(t2) = true;
            
            tb1 = find_neib(tc, t1);
            tb2 = find_neib(tc, t2);
            tb = union(tb1, tb2);
            % disconnect neighboring triangles
            tc(t1, :) = 0; tc(:, t1) = 0;
            tc(t2, :) = 0; tc(:, t2) = 0;
            % reconnect external triangles
            for tn = tb
                if is_neighbour(triangles, tn, t1)
                    tc = connecttril(tc, tn, t1);
                    triangles_changed(tn) = 1;
                end
                if is_neighbour(triangles, tn, t2)
                    tc = connecttril(tc, tn, t2);
                    triangles_changed(tn) = 1;
                end
            end
            tc = connecttril(tc, t1, t2);
        end
    end
end

function tsd = tsetdiff(triangles, t1, t2)
    tsd = find(triangles(:,t1) & ~triangles(:, t2));
end

function tc = connecttril(tc, u,v)
  tc(max(u, v), min(u, v)) = 1;
end

function ti = tintersect(triangles, t1, t2)
    ti = find(triangles(:, t1) & triangles(:, t2));
end

function ti = txor(triangles, t1, t2)
    ti = find(xor(triangles(:, t1),triangles(:, t2)));
end

function tn = find_neib(tc, t1)
  tn = union(find(tc(t1,:)), find(tc(:,t1)));
  if size(tn, 1) ~= 1, tn = tn'; end
end