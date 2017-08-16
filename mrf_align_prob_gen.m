function [Mask] = mrf_align_prob_gen(G_s, G_t, Cost, X_hub_align)
%
lambda = 0.01;
[numS, numT] = size(Cost);
X = ones(numS, numT)/sqrt(numT);


for iter = 1:20
    X = Cost+ 1e-5 + G_s*X*G_t*lambda;
    s = 1./sqrt(sum(X'.*X'))';
    X = (s*ones(1,numT)).*X;
end

Mask = X_hub_align;
Mask = max(Mask, double(X > 1e-2));
Mask = max(Mask, Cost);