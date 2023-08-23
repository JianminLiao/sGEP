function [x, f] = TRF_modified(A,B,eta,options, x0)
%% Set Parameters
cardinality=options.cardinality;
optTol=options.optTol;
maxIter=options.maxIter;
%% Default initialization
% if nargin < 5
%     [~,idx]=sort(diag(A), 'descend');
%     x0 = zeros(size(A,1),1);
%     x0(idx(1:cardinality)) = 1;
%     x0 = x0 / norm(x0);
% end

% if nargin < 5
%     [x0,~]=eigs(A,B,1,'largestreal');
% end

% if nargin < 5
%     x0 = randn(length(A),1);    
%     x0  = x0/norm(x0);
% end
% global VALUE
x = truncate_operator(x0, cardinality);
n=length(A);
% power step
Ax=A*x;
Bx=B*x;
rho=(x'*Ax)/(x'*Bx);
I=eye(n);
C=I+eta/rho*A-eta*B;
s = C*x;
f = rho;

% truncate step
x = truncate_operator(s, cardinality);
f_old = f;
rho_old=rho;

i=1;
%% Main algorithmic loop
while i<=maxIter
    % power step
    Ax=A*x;
    Bx=B*x;
    rho=(x'*Ax)/(x'*Bx);
    C=C+eta*(1/rho-1/rho_old)*A; %faster
    %C=I+(eta/rho)*(A-rho*B);
    s = C*x;
    
    % truncate step
    x = truncate_operator(s, cardinality);
    f = rho;
%     VALUE=[VALUE;f];
    if (f - f_old)/f_old < optTol
        break;
    end
    f_old = f;
    rho_old=rho;
    i=i+1;
end
end

%% Evaluate the truncate operator
function xc = truncate_operator(x,r)
[~,idx]=maxk(x,r,'ComparisonMethod','abs');
xc=sparse(length(x), 1);
x_restrict = x(idx);
xc(idx) = x_restrict / norm(x_restrict);
end
