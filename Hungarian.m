function [rowsol] = Hungarian(Perf)

[ns, nt] = size(Perf);
Perf = (Perf - min(min(Perf)))/(max(max(Perf)) - min(min(Perf)));
if ns <= nt
    CostMat = ones(nt, nt)*16384;
    CostMat(1:ns, :) = floor(Perf*8192);
    rowsol = lapjv(CostMat, 0);
    rowsol = rowsol(1:ns);
else
    CostMat = ones(ns, ns)*16384;
    CostMat(1:nt,:) = floor(Perf*8192)';
    rowsol = lapjv(CostMat, 0);
    rowsol = rowsol(1:nt);
end

