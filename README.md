# This code implements a fast graph matching problem, which relaxes quadratic assignment.
oad('test_data.mat');
load('parameters.mat');
% Without edge similarity matrix
[Corres] = graph_matching_master(Graph_s, Graph_t, NodeSimilarity, Para);
% With edge similarity matrix
[Corres] = graph_matching_master2(edges_s, edges_t, NodeSimilarity, EdgeSimilarity, Para);
