function [Corres] = graph_matching_master(Graph_s, Graph_t, NodeSimilarity, Para)
% This function returns an injective map from graph_s to graph_t. We assume
% that the number of vertices of the source graph is smaller than that of
% the target graph.
% Input arguments:
%   Graph_s: The source graph, which is given by a sparse matrix
%   Graph_t: The target graph, which is also given by a sparse matrix
%   NodeSimilarity: A sparse matrix the specifies the vertex similarities
%                 between vertices of the source graph and vertices of the
%                 target graph
%   Para: parameters
%         Para.lambda_edge: strength of the pair-wise term
%         Para.lambda_mu: a large constant enforcing the consistentcy of
%                         bi-directional correspondences
% Output:
%   Corres: a matrix of dimension 2x|V_s|, specifiying vertex indices of
%           correspondences
if isfield(Para, 'flag_fast') == 1 && Para.flag_fast == 1
    X_st = mrf_align_admm2(Graph_s, Graph_t, NodeSimilarity,...
        Para.lambda_edge,...
        Para.mu);
else
    Data = mrf_align_opt_data(Graph_s, Graph_t, NodeSimilarity);
    X_st = mrf_align_admm(Data);
end

% Round the fractional solution into an integer solution
ns = size(Graph_s, 1);
nt = size(Graph_t, 1);
rowsol = Hungarian(1-X_st);
Corres = [1:ns; rowsol];