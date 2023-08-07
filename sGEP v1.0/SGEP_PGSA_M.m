function [x, G, iter] = SGEP_PGSA_M(A,B,r,reltol,alpha0,x)

%% The code is associated with the article "FIRST-ORDER ALGORITHMS FOR A CLASS OF FRACTIONAL OPTIMIZATION PROBLEM"
%% PGSA_L(N=0) for SGEP 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       min  G(x) := (x^T B x) / (x^T A x)                            			%%%
%%%		  s.t.  |x|_0  <=  r  														%%%
%%%                                                                     			%%%
%%% Input: 	A and B: matrices of the problem			            				%%%
%%%			r: a positive integer, sparsity constraint								%%%
%%%			x: initial point, which satisfies |x|_0 <= r and  (x^T A x) neq 0		%%%
%%%			alpha0: the minimum step size (alpha0 = 0.99/norm(B,2) is recommended)	%%%
%%%			reltol: relative tolerance     											%%%
%%%			maxit:  the upper limit of iteration number								%%%
%%%																					%%%
%%% Output: x: the last iteration point	xn											%%%
%%%			G: the value of G(xn)                                           		%%%
%%%			iter: iteration number in fact											%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% parameters
maxit = 6000;	% max iteration times
a = 1e-3;		% the same as 'a' in Alg2 
rho = 0.5;		% the same as '\eta' in Alg2
alpha = alpha0;		% 'alpha' in the first iteration


%% iteration
iter = 0;
x=sparse(x);
Bx = B*x;
Ax = A*x;
G = x'*Bx/(x'*Ax);
stop_criterion = false;

while(~stop_criterion && iter < maxit)
	iter = iter +1;
	D = -Bx+G*Ax;
	G_old = G;
    x_old = x;
    Bx_old = Bx;	
	
	% line-search process
	line_search_stop_condition = false;	
	while ~line_search_stop_condition
		x = ProjC(x_old+alpha*D,r);
        Bx = B*x;
        Ax = A*x;
        G = x'*Bx/(x'*Ax);
		diff_x = x-x_old;
		diff_x_norm22 = diff_x'*diff_x;
		line_search_stop_condition = (alpha==alpha0 | G <= G_old - a/2*diff_x_norm22);
		alpha = max(alpha*rho,alpha0);
	end
	alpha = max(alpha0,diff_x_norm22 / (diff_x' * (Bx-Bx_old)));
	
	stop_criterion = norm(x-x_old,2)/norm(x,2) < reltol;
end
	
end

function xc = ProjC(x,r)
[~,idx]=maxk(x,r,'ComparisonMethod','abs');
xc=sparse(length(x), 1);
x_restrict = x(idx);
xc(idx) = x_restrict / norm(x_restrict);
end