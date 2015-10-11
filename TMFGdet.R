#function [A,tri,separators,cliques,cliqueTree]=TMFGdet(W)

detgain <- function(vid, tid, W, tri){
  k4 <- simplify2array(c(tri[tid,], vid))
  k3 <- simplify2array(c(tri[tid,]))
  det(W[k3,k3]) / det(W[k4,k4])
}

mdetgain <- function(ous, tris, W, tri){
  mapply(detgain, ous, tris, MoreArgs=list(W, tri))
}


TMFGdet <- function(W){
  # N    = size(W,1);
  # A    = sparse(N,N);     % ininzialize adjacency matrix
  # in_v = zeros(N,1);      % ininzialize list of inserted vertices
  # tri  = zeros(2*N-4,3);  % ininzialize list of triangles
  # separators=zeros(N-4,3);% ininzialize list of 3-cliques (non face-triangles)
  
  N <- nrow(W)
  A <- rep(0, N * N)
  dim(A) <- c(N, N)
  in_v <- rep(0, N)
  tri <- rep(0, 3*(2*N-4))
  dim(tri) <- c(2*N-4,3)
  separators <- rep(0,(N-4)*3)
  dim(separators) <- c(N-4,3)
# %% find 3 vertices with largest strength
# s    = sum(W.*(W>mean(W(:))),2);
# [~,j]=sort(s,'descend');
# in_v(1:4)  = j(1:4);
# ou_v = setdiff([1:N],in_v); % list of vertices not inserted yet
  X <- W * (W > mean(W) )
  s <- rowSums(X)
  k4 <- order(s, decreasing=T)[1:4]
  in_v[1:4] <- k4
  ou_v <- setdiff(1:N, k4)

# %% build the tetrahedron with largest strength
# tri(1,:)=in_v([1 2 3]);
# tri(2,:)=in_v([2 3 4]);
# tri(3,:)=in_v([1 2 4]);
# tri(4,:)=in_v([1 3 4]);
  tri[1,] <- in_v[c(1,2,3)]
  tri[2,] <- in_v[c(2,3,4)]
  tri[3,] <- in_v[c(1,2,4)]
  tri[4,] <- in_v[c(1,3,4)]
# A(in_v(1),in_v(2)) = 1; 
# A(in_v(1),in_v(3)) = 1;
# A(in_v(1),in_v(4)) = 1;
# A(in_v(2),in_v(3)) = 1;
# A(in_v(2),in_v(4)) = 1;
# A(in_v(3),in_v(4)) = 1;
  A[in_v, in_v] <- 1; diag(A) <- 0;
# %% build initial gain table
# gain = -inf(N,2*N-4);
# for vect=ou_v
# gain(vect,1) = (det(W(tri(1,:),tri(1,:)))/det(W([vect tri(1,:)],[vect tri(1,:)])));
# gain(vect,2) = (det(W(tri(2,:),tri(2,:)))/det(W([vect tri(2,:)],[vect tri(2,:)])));
# gain(vect,3) = (det(W(tri(3,:),tri(3,:)))/det(W([vect tri(3,:)],[vect tri(3,:)])));
# gain(vect,4) = (det(W(tri(4,:),tri(4,:)))/det(W([vect tri(4,:)],[vect tri(4,:)])));
# %     gain(vect,1) = sum(W(vect,tri(1,:)).^2,2);% approx (works better!)
# %     gain(vect,2) = sum(W(vect,tri(2,:)).^2,2);% approx (works better!)
# %     gain(vect,3) = sum(W(vect,tri(3,:)).^2,2);% approx (works better!)
# %     gain(vect,4) = sum(W(vect,tri(4,:)).^2,2);% approx (works better!)
# end
  gain <- rep(-Inf, N* (2*N-4)); dim(gain) <- c(N, 2*N-4)
  gain[ou_v, 1:4] <- outer(ou_v, 1:4, mdetgain,  W, tri)
  kk <- 4
  maxval <- apply(gain, 2, max)
  bestv  <- apply(gain, 2, which.max)
# kk = 4;  % number of triangles
# for k=5:N
  for (k in 5:N) {
    # %% find best vertex to add in a triangle
    # if length(ou_v)==1 %special case for the last vertex
    # ve = ou_v;
    # v  = 1;
    # [~,tr] = max(gain(ou_v,:));
    # else
    #   [gij,v]= max(gain(ou_v,:));
    # [~,tr] = max( gij );
    # ve = ou_v(v(tr));
    # v  = v(tr);
    # end
    tr <- which.max(maxval)
    ve <- bestv[t]
    # %% update vertex lists
    # ou_v = ou_v([1:(v-1),(v+1):end]);
    # in_v(k)   = ve;
    ou_v <- ou_v[ou_v != v]
    in_v[k] <- ve
    # %% update adjacency matrix
    # A(ve,tri(tr,:))=1;
    A[ve, tri[tr,]] <- 1; A[tri[tr,], ve] <-1
    # %% update 3-clique list
    # separators(k-4,:) = tri(tr,:);
    separators[k-4, ] <- tri[tr, ]
    # %% update triangle list replacing 1 and adding 2 triangles 
    # tri(kk+1,:) = [tri(tr,[1,3]),ve]; % add
    # tri(kk+2,:) = [tri(tr,[2,3]),ve]; % add
    # tri(tr,:)   = [tri(tr,[1,2]),ve]; % replace
    tri[kk+1,] <- simplify2array(c(tri[tr,1], tri[tr,3], ve))
    tri[kk+2,] <- simplify2array(c(tri[tr,2], tri[tr,3], ve))
    tri[tr  ,] <- simplify2array(c(tri[tr,1], tri[tr,2], ve))
    # %% update gain table
    # gain(ve,:)=0;
    # for vect=ou_v
    # gain(vect,tr)  = (det(W(tri(tr,:)  ,tri(tr,:)))  /det(W([vect,tri(tr,:)]  ,[vect,tri(tr,:)])))  ;
    # gain(vect,kk+1)= (det(W(tri(kk+1,:),tri(kk+1,:)))/det(W([vect,tri(kk+1,:)],[vect,tri(kk+1,:)]))); 
    # gain(vect,kk+2)= (det(W(tri(kk+2,:),tri(kk+2,:)))/det(W([vect,tri(kk+2,:)],[vect,tri(kk+2,:)]))); 
    # %       gain(vect,tr)  = sum(W(vect,tri(tr,:)).^2,2)  ;% approx (works better!)
    # %     	gain(vect,kk+1)= sum(W(vect,tri(kk+1,:)).^2,2);% approx (works better!)
    # %       gain(vect,kk+2)= sum(W(vect,tri(kk+2,:)).^2,2);% approx (works better!) 
    # end;
    gain[ve, ] <- -Inf
    newt <- simplify2array(c(tr, kk+1, kk+2 ))
    gain[ou_v, newt] <- outer(ou_v, newt, mdetgain,  W, tri)
    # %% update number of triangles
    # kk = kk+2; 
    # if mod(k,1000)==0,fprintf('TMFG T2 det: %0.2f per-cent done\n',k/N*100);end
    # end
    kk <- kk + 2 
  }
  A <- W * as.vector (A > 0)
  # %% comutes 4-clique list
  # if nargout>3     
  #     cliques = [in_v(1:4)';separators,in_v(5:end)]; 
  # end
  cliques <- rbind(t(in_v[1:4]),
                   cbind(separators, in_v[5:N]))
  A
}
