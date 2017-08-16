function [X] = mrf_align_admm(Data)
% This function implements an ADMM method for solving the linear
% programming relaxation of the MRF problem described in Data
% The linear programming relaxation is described below (the dual variables
%                                                       are on the right-hand side)
% maximize     <W_X, X > + <w_y, y>
% subject to   y     >= 0                        : z0  >= 0
%              E_i*y <= s                        : z1i >= 0
%              s     == F_i*X                    : z2i >= 0
%              X* 1  <= 1                        : z3  >= 0
%              X'*1  <= 1                        : z4  >= 0
%              X     >= 0                        : Z5  >= 0
%
% the augmented Lagrangian to be minimized is given by
% <1, z3> + <1, z4> + <z1 - z2, s> + <w_y + z0 - E'*z1, y> 
% + <W_X + Z5 + F'*z2 - z3*1'-1*z4', X> 
% + 1/2u*(\|z1 - z2\|^2 + \|w_y + z0 - E'*z1\|^2 + \|W_X + Z5 + F'*z2 - z3*1'-1*z4'\|^2)
% 

% Making sure that the cost function is non-degenerate
[numV_s, numV_t] = size(Data.W_X);
dimX = numV_s*numV_t;
dimy = length(Data.w_y);

X = zeros(numV_s, numV_t); % Store the solution matrix
y = zeros(dimy, 1); % Store the y vector
s = zeros(length(Data.idsF), 1); % primal variable s

z0 = y; % dual variable of constraint y >= 0
z1 = s; % dual variable of constraint Ey <= s
z2 = z1; % dual variable of constraint s = \set{F}(X)
z3 = zeros(numV_s, 1); % dual variable of constraint X1 <= 1
z4 = zeros(numV_t, 1); % dual variable of constraint X'1 <= 1
Z5 = zeros(numV_s, numV_t); % dual variable of constraint X >= 0

% pre-compute the number of non-zeros elements in each row of E
dE = sum(Data.E')' + 1;
% group the rows of E so that in each sub-matrix, there is only one
% non-zero entry per column
ids1 = (Data.offs(1)+1):Data.offs(2);
ids2 = (Data.offs(2)+1):Data.offs(3);
ids3 = (Data.offs(3)+1):Data.offs(4);
ids4 = (Data.offs(4)+1):Data.offs(5);
% Buffer the matrices
E1 = Data.E(ids1, :);
E2 = Data.E(ids2, :);
E3 = Data.E(ids3, :);
E4 = Data.E(ids4, :);
%
Et = Data.E';
E1t = Data.E(ids1, :)';
E2t = Data.E(ids2, :)';
E3t = Data.E(ids3, :)';
E4t = Data.E(ids4, :)';
%
dE1 = dE(ids1);
dE2 = dE(ids2);
dE3 = dE(ids3);
dE4 = dE(ids4);
%

% find the elements of X (e.g., the correspondences) that are really being
% constrainted by the pairwise constraints
numCons = length(Data.idsF);
F = sparse(1:numCons,...
    Data.idsF,...
    ones(1, numCons),...
    numCons, dimX);
Ft_full = F';

constrainedIds = find(sum(F) > 0);
F = F(:, constrainedIds);
% Buffer matrices
Ft = F';
dF = (Ft*ones(size(F,1),1)+1);
%
mu = 1;
rho = 1.001;

fprintf(' [');
for iteration = 1:4000
    % Optimize z0
    z0 = max(0, Et*z1 - mu*y - Data.w_y);
    
    % Optimize z1
    tp1 = Data.w_y + z0 + mu*y;
    tp2 = z2 - mu*s;
    %
    tp21 = tp2(ids1);
    tp22 = tp2(ids2);
    tp23 = tp2(ids3);
    tp24 = tp2(ids4);
    %
    for innerIter = 1:1
        % Optimize z11
        TP1 = tp1 - Et*z1 + E1t*z1(ids1);
        z1(ids1) = max(0, (E1*TP1 + tp21)./dE1);
    
        % Optimize z12 
        TP1 = tp1 - Et*z1 + E2t*z1(ids2);
        z1(ids2) = max(0, (E2*TP1 + tp22)./dE2);
        
        % Optimize z13
        TP1 = tp1 - Et*z1 + E3t*z1(ids3);
        z1(ids3) = max(0, (E3*TP1 + tp23)./dE3);
        
        % Optimize z14 
        TP1 = tp1 - Et*z1 + E4t*z1(ids4);
        z1(ids4) = max(0, (E4*TP1 + tp24)./dE4);
    end
    
    % Optimize z2
    TP1 = z3*ones(1,numV_t) + ones(numV_s,1)*z4' - Data.W_X -mu*X - Z5;
    TP1 = TP1(constrainedIds)';
    tp1 = z1 + mu*s;
    r = (Ft*tp1 - TP1)./dF;
    z2 = tp1 - F*r;
    
    % Optimize the z3 and z4 in order
    TP2 = Data.W_X + reshape(full(Ft_full*z2), [numV_s, numV_t]); 
    TP1 = TP2 + Z5 + mu*X;
    
    tp1 = TP1*ones(numV_t,1);
    tp2 = TP1'*ones(numV_s, 1);
    for inner = 1:5
        z3 = max(0, (tp1 - (sum(z4) + mu))/numV_t);
        z4 = max(0, (tp2 - (sum(z3) + mu))/numV_s);
    end
    
    % Optimize Z5
    TP1 = z3*ones(1, numV_t) + ones(numV_s,1)*z4' - TP2;
    Z5 = max(0, TP1 - mu*X);
    
    % Update the primal variable X
    X = X + (Z5-TP1)/mu;
    
    % Update the primal variable y
    y = y + (Data.w_y + z0 - Et*z1)/mu;
    
    % Update the primal variable s
    s = s + (z1-z2)/mu;
    
    mu = mu*rho;
    if mod(iteration, 40) == 0
        fprintf(' iter = %d, [%d,%d], energy = [%.3f, %.2f], norm(X) = [%.3f, %.3f]\n', iteration, numV_s, numV_t, sum(sum(Data.W_X.*X)),  Data.w_y'*y,  max(sum(X')), min(sum(X')));
    end
end
fprintf(']\n');
