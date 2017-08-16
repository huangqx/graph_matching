function [X_st] = mrf_align_admm2(G_s, G_t, W_st, lambda, mu)
% A faster and simpler admm implementation
% We do appropriate normalizations
% W_st = normalize_score(W_st);
G_s = normalize_graph(G_s);
G_t = normalize_graph(G_t);

%
[ns, nt] = size(W_st);
%
C_st = zeros(ns, nt+1);
C_st(:,1:nt) = full(W_st);
%
C_ts = zeros(nt, ns+1);
C_ts(:, 1:ns) = full(W_st)';
%
X_st = ones(ns, nt+1)/(nt+1);
X_ts = ones(nt, ns+1)/(ns+1);
%
for iter = 1:30
    TP1 = lambda*G_s*X_st(:,1:nt)*G_t + mu*X_ts(:,1:ns)'; 
    TP2 = lambda*G_t*X_ts(:,1:ns)*G_s + mu*X_st(:,1:nt)';
    
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
    TP1 = lambda*G_s*X_st(:,1:nt)*G_t + mu*X_ts(:,1:ns)'; 
    TP2 = lambda*G_t*X_ts(:,1:ns)*G_s + mu*X_st(:,1:nt)';
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
X_st = X_st(:, 1:nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize the similarity score between the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize the similarity graph among the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G_nor] = normalize_graph(G)
% We will high-weights to edges 
d = max(full(sum(G)), 1);
ids_non_hub = find(d < max(d)/12);
[rows, cols, vals]= find(G);
d = d/max(d);
vals_nor = max(d(rows'), d(cols'));
if 1
    vals_nor = sqrt(d(rows').*d(cols'));
end

n = size(G, 1);
G_nor = sparse(rows, cols, vals_nor, n, n);
G_nor(ids_non_hub, ids_non_hub) = 0;
G_nor = sparse(G_nor);