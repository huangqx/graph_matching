# This code implements a fast graph matching problem, which relaxes quadratic assignment.
load('test_data.mat');

load('parameters.mat');

% Without edge similarity matrix (edge similarity scores are 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Corres = graph_matching_master(Graph_s, Graph_t, NodeSimilarity, Para);

% With edge similarity matrix (edge similarity scores are varying)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Corres = graph_matching_master2(edges_s, edges_t, NodeSimilarity, EdgeSimilarity, Para);
