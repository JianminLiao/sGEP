function [x, f] = TPower_modified(A, options, x0)
%%  Truncated power method for sparse eigenvalue problem.
%
%     max x'*A*x    subject to ||x||=1, ||x||_0 <= k.
%
% *** inputs ***
% - A:               p x p symmetric positive semidefinite matrix
% - options:         a structure stores user-specified parameters which include:
%    -- cardinality: cadinality k
%    -- optTol:      optimality tolerance (default: 1e-6)
%    -- maxIter:     max number of iteration (default: 50)
% - x0:              initialization vector
%
% *** outputs ***
% - x:            p-dimensional sparse eigenvector vector with k non-zeros
% - f:            objective value at the output x
%
%   Modified from
%   Xiao-Tong Yuan, Tong Zhang, Truncated Power Method for Sparse Eigenvalue Problems, Technical Report, 2011
%
% Copyright (C) 2011/2012 - Xiao-Tong Yuan.

%% Set Parameters
cardinality=options.cardinality;
optTol=options.optTol;
maxIter=options.maxIter;

%% Default initialization
if nargin < 3
    [~,idx]=sort(diag(A), 'descend');
    x0 = zeros(size(A,1),1);
    x0(idx(1:cardinality)) = 1;
    x0 = x0 / norm(x0);
end

x = sparse(x0);
x = truncate_operator(x, cardinality);
% power step
s = A*x;
g = 2*s;
f = x'*s;

% truncate step
x = truncate_operator(g, cardinality);
f_old = f;

i=1;
%% Main algorithmic loop
while i<=maxIter
    % power step
    s = A*x;
    g = 2*s;
    
    % truncate step
    x = truncate_operator(g, cardinality);
    f = x'*s;
    
    if (f - f_old)/f_old < optTol
        break;
    end
    f_old = f;
    i=i+1;
end
end

%% Evaluate the truncate operator
function u = truncate_operator(v , k)
u = zeros(length(v), 1);
[~,idx]=maxk(abs(v),k);
v_restrict = v(idx);
u(idx(1:k)) = v_restrict / norm(v_restrict);
u = sparse(u);
end
