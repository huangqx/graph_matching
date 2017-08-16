# This code implements a fast graph matching problem, which relaxes quadratic assignment.
>>load('test_data.mat');

>>load('parameters.mat');

>>[Corres] = graph_matching_master(Graph_s, Graph_t, NodeSimilarity, Para);
