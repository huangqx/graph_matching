function [Corres] = graph_matching_master2(edges_s,...
    edges_t,...
    NodeSimilarity,...
    EdgeSimilarity,...
    Para)
% This function returns an injective map from graph_s to graph_t. We assume
% that the number of vertices of the source graph is smaller than that of
% the target graph.
% Input arguments:
%   edges_s: edges of the source graph, which is given by a matrix of 
%            dimension 2 x numEdges_sourceGraph
%   edges_t: edges of the target graph, which is given by a matrix of 
%            dimension 2 x numEdges_targetGraph
%   NodeSimilarity: A sparse matrix that specifies the vertex similarities
%                 between vertices of the source graph and vertices of the
%                 target graph
%   EdgeSimilarity: A sparse matrix that specifies the edge similarities 
%                 between edges of the source graph and edges of the target
%                 graph
%   Para: parameters
%         Para.lambda_edge: strength of the pair-wise term
%         Para.lambda_mu: a large constant enforcing the consistentcy of
%                         bi-directional correspondences
% Output:
%   Corres: a matrix of dimension 2x|V_s|, specifiying vertex indices of
%           correspondences
if isfield(Para, 'flag_fast') == 1 && Para.flag_fast == 1
    X_st = mrf_align_admm3(edges_s, edges_t, NodeSimilarity, EdgeSimilarity,...
        Para.lambda_edge,...
        Para.mu);
else
    Data = mrf_align_opt_data2(edges_s, edges_t, NodeSimilarity, EdgeSimilarity);
    X_st = mrf_align_admm(Data);
end

% Round the fractional solution into an integer solution
[ns, nt] = size(NodeSimilarity);
rowsol = Hungarian(1-X_st);
Corres = [1:ns; rowsol];