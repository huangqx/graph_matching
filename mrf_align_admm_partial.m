function [X_st] = mrf_align_admm_partial(edges_s, edges_t,...
    W_node_st, W_edge_st,...
    lambda,...
    mu,...
    gamma)
% A faster and simpler admm implementation
% We do appropriate normalizations
% lambda: tradeoff parameter for the edge potential
% mu:     tradeoff parameter for enforcing bidirectional consistency 
% gamma:  a hyper-parameter used to remove correspondences that are weak,
%         i.e., node similarity is low, and 
% Restore the pair-wise term in the matrix-form (a sparse matrix)
W_pair = pairwise_term_generator(edges_s, edges_t, W_edge_st);
%
[ns, nt] = size(W_node_st);
%
C_st = zeros(ns, nt+1);
C_st(:,1:nt) = full(W_node_st);
C_st(:, nt+1) = gamma + mu;
%
C_ts = zeros(nt, ns+1);
C_ts(:, 1:ns) = full(W_node_st)';
C_ts(:, ns+1) = gamma + mu;
%
X_st = ones(ns, nt+1)/(nt+1);
X_ts = ones(nt, ns+1)/(ns+1);
%
for iter = 1:30
    TP1 = mu*X_ts(:,1:ns)';
    TP1 = TP1 + lambda*pairwise_term_multiply(W_pair, X_st(:,1:nt));
    TP2 = mu*X_st(:,1:nt)';
    TP2 = TP2 + lambda*pairwise_term_multiply(W_pair, X_ts(:,1:ns)')';
    %
    X_st = C_st;
    X_st(:, 1:nt) = X_st(:, 1:nt) + TP1;
    X_ts = C_ts;
    X_ts(:, 1:ns) = X_ts(:, 1:ns) + TP2;
    
    % Normalize
    s = sqrt(sum(X_st'.*X_st')');
    X_st = X_st./kron(s, ones(1,nt+1));
    
    s = sqrt(sum(X_ts'.*X_ts')');
    X_ts = X_ts./kron(s, ones(1,ns+1));
end
%
beta_min = 0.001;
beta_max = 0.01;
t = (0:29)/29;

for iter = 1:30
    beta = exp(log(beta_min)*(1-t(iter))+log(beta_max)*t(iter));
    TP1 = mu*X_ts(:,1:ns)'; 
    TP1 = TP1 + lambda*pairwise_term_multiply(W_pair, X_st(:,1:nt));
    TP2 = mu*X_st(:,1:nt)';
    TP2 = TP2 + lambda*pairwise_term_multiply(W_pair, X_ts(:,1:ns)')';
    %
    X_st_2 = C_st;
    X_st_2(:, 1:nt) = X_st_2(:, 1:nt) + TP1;
    X_ts_2 = C_ts;
    X_ts_2(:, 1:ns) = X_ts_2(:, 1:ns) + TP2;
    %
    X_st_2 = exp(X_st_2*beta);
    X_ts_2 = exp(X_ts_2*beta);
    X_st = X_st.*X_st_2;
    X_ts = X_ts.*X_ts_2;
    % Normalize
    s = sum(X_st')';
    X_st = X_st./kron(s, ones(1,nt+1));
    
    s = sum(X_ts')';
    X_ts = X_ts./kron(s, ones(1,ns+1));
end
% Determine the vertices that do not have correspondences
rowBounds = X_st(:, nt+1);
X_st = X_st(:, 1:nt);
ids = find(max(X_st') < rowBounds');
X_st(ids, :) = 0; 

% Generate the pair-wise term
function [W_pair] = pairwise_term_generator(edges_s, edges_t, W_edge_st)
%
nv_s = max(max(edges_s));
nv_t = max(max(edges_t));
dim = nv_s*nv_t;
[edgeCorres_s_ids, edgeCorres_t_ids, edgeCorres_weights] = find(W_edge_st);
% We reuse variables to save space
edgeCorres_s_ids = edges_s(:, edgeCorres_s_ids); % vertices of the source edge
edgeCorres_t_ids = edges_t(:, edgeCorres_t_ids); % vertices of the target edge
% indices of the source correspondence
edgeCorres_s_ids(1,:) = (edgeCorres_t_ids(1,:)-1)*nv_s + edgeCorres_s_ids(1,:);
% indices of the target correspondence
edgeCorres_s_ids(2,:) = (edgeCorres_t_ids(2,:)-1)*nv_s + edgeCorres_s_ids(2,:);
W_pair = sparse([edgeCorres_s_ids(1,:);edgeCorres_s_ids(2,:)],...
    [edgeCorres_s_ids(2,:);edgeCorres_s_ids(1,:)],...
    kron(ones(2,1), edgeCorres_weights),...
    dim, dim);

%
function [X_st_out] = pairwise_term_multiply(W_pair, X_st)
%
[ns, nt] = size(X_st);
X_st = reshape(X_st, [ns*nt,1]);
X_st_out = W_pair*X_st;
X_st_out = reshape(X_st_out, [ns, nt]);