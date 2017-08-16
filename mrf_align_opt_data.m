function [Data] = mrf_align_opt_data(G_s, G_t, Blast)
% Generate the linear programming data for graph matching
% Input arguments:
%       G_s: the source input graph
%       G_t: the target input graph
%     Blast: the computed pair-wise blast scores
G_s = max(G_s, G_s');
G_t = max(G_t, G_t');
%
d_s = full(sum(G_s));
d_t = full(sum(G_t));
%
ids_nonhub_s = find(d_s < 6);
ids_nonhub_t = find(d_t < 6);
G_s(ids_nonhub_s, ids_nonhub_s) = 0;
G_t(ids_nonhub_t, ids_nonhub_t) = 0;
%
[rows_S, cols_S, vals_S] = find(G_s);
edges_S = [cols_S'; rows_S'];
%
[rows_T, cols_T, vals_T] = find(G_t);
edges_T = [cols_T'; rows_T'];
%
ns = size(G_s, 1);
nt = size(G_t, 1);
[rows, cols, vals] = find(ones(ns, nt));
mask_data = [rows';cols'];
[E1_rowIds, E2_rowIds, rowCorIds1] = BinaryCons(edges_S, edges_T, mask_data, [ns,nt]);
mask_data = [cols';rows'];
[E3_rowIds, E4_rowIds, rowCorIds2] = BinaryCons(edges_T, edges_S, mask_data, [nt,ns]);

dimCor = length(rowCorIds1);
dimY = size(E1_rowIds, 2);
E1 = sparse(E1_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF1 = find(sum(E1')>0);
E1 = E1(idsF1, :);
idsF1 = rowCorIds1(idsF1);

E2 = sparse(E2_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF2 = find(sum(E2')>0);
E2 = E2(idsF2, :);
idsF2 = rowCorIds1(idsF2);

dimCor = length(rowCorIds2);
dimY = size(E3_rowIds, 2);
E3 = sparse(E3_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF3 = find(sum(E3')>0);
E3 = E3(idsF3, :);
idsF3 = rowCorIds2(idsF3);

E4 = sparse(E4_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF4 = find(sum(E4')>0);
E4 = E4(idsF4, :);
idsF4 = rowCorIds2(idsF4);

TP = reshape(1:(ns*nt), [ns, nt])';
TP = reshape(TP, [1, ns*nt]);
idsF3 = TP(idsF3);
idsF4 = TP(idsF4);

Mask = reshape(1:(ns*nt), [ns, nt]);

[v1Ids, s, v] = find(E1);
[v2Ids, s, v] = find(E2);
[v3Ids, s, v] = find(E3);
[v4Ids, s, v] = find(E4);
v1Ids = Mask(idsF1(v1Ids'));
v2Ids = Mask(idsF2(v2Ids'));
v3Ids = Mask(idsF3(v3Ids'));
v4Ids = Mask(idsF4(v4Ids'));
num = length(v1Ids);

M1 = sparse(min(v1Ids, v2Ids), max(v1Ids, v2Ids), 1:num);
M2 = sparse(min(v3Ids, v4Ids), max(v3Ids, v4Ids), 1:num);
[v1Ids, v2Ids, order12] = find(M1);
[v3Ids, v4Ids, order34] = find(M2);
E1 = E1(:, order12');
E2 = E2(:, order12');
E3 = E3(:, order34');
E4 = E4(:, order34');


Data.E = [E1;E2;E3;E4];
Data.offs = [0, length(idsF1), length(idsF2), length(idsF3), length(idsF4)];
Data.idsF = [idsF1, idsF2, idsF3, idsF4];
for i = 1:4
    Data.offs(i+1) = Data.offs(i) + Data.offs(i+1);
end
Data.w_y = ones(size(Data.E, 2), 1);
Data.W_X = Blast;
% 
% Data.dimX = size(Cost, 1)*size(Cost, 2);
% Data.dimY = dimY;
% 
% Data.G_s = G_s;
% Data.G_t = G_t;

% len = length(Data.idsF);
% Data.Ft = sparse(Data.idsF, 1:len, ones(1, len),...
%     Data.dimX, len);
%[s, idss] = sort(-sum(double(G_s)));
%[s, idst] = sort(-sum(double(G_t)));

% score_s = score_s/mean(score_s);
% score_t = score_t/mean(score_t);
% [ns, nt] = size(Cost);
% ds = score_s'*ones(1, nt);
% dt = ones(ns,1)*score_t;

% Data.D_Cost = min(ds, dt)./power(max(ds,dt), 0.33);
% Data.D_Cost = Data.D_Cost/max(max(Data.D_Cost));

