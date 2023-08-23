function [x, G, iter] = SGEP_IFTRR(A,B,r,x,reltol,maxit,norm_A,norm_B)
%% The code is associated with the article "An Inverse-free Truncated Rayleigh-Ritz Method for Sparse Generalized Eigenvalue Problem"
%% Authors of the paper: Yunfeng Cai and Ping Li
%% IFTRR for SGEP 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       max  rho(x) := (x^T A x) / (x^T B x)                            			%%%
%%%		  s.t.  |x|_0  <=  r  														%%%
%%%                                                                     			%%%
%%% Input: 	A and B: matrices of the problem			            				%%%
%%%			r: a positive integer, sparsity constraint								%%%
%%%			x: initial point, which satisfies |x|_0 <= r 							%%%
%%%			reltol: |ρt - ρt-1| < reltol                                          %%%
%%%			maxit:  the upper limit of outer iteration number						%%%
%%%																					%%%
%%% Output: x: the last iteration point	xn											%%%
%%%			G: the value of 1/rho(xn)                                          		%%%
%%%			outer_iter: outer iteration number in fact								%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% parameters
delta_k = 30;
m = 5;
tol  = 1e-3; 
tol1 = 0.01;
tol3 = 1e-9; 

s2 = r+delta_k;
Dimension = length(x);
iter = 0;
rho_iter = x'*A*x/ (x'*B*x);
rho_iter_old=rho_iter;
stop_criterion = false;

while(~stop_criterion && iter < maxit)
	iter = iter+1;
	[Q,~] = Arnoldi0(A-rho_iter*B,m-1,x);
    [tilde_y,~] = leadingEig(Q'*A*Q,Q'*B*Q);
	w = Q*tilde_y;
	[~,I]=sort(abs(w),'descend');
    
    J = false(Dimension,1);
    J(I(1:s2)) = true;
    J = ill_condition_process(B,J,tol3);
	[z,rho_s2]=leadingEig(A(J,J),B(J,J));
    
    rho_iter = rho_s2;
    a=r;b=s2;
    while b>a
        s = a+floor((b-a)/2);
        J = false(Dimension,1);
        J(I(1:s)) = true;  
        J = ill_condition_process(B,J,tol3);
        [z,rho_iter]=leadingEig(A(J,J),B(J,J));
        if rho_s2-rho_iter<=tol
            b=s;
        else
            a=s+1;
        end
    end
	
	
	x_old = x;
	x = zeros(Dimension,1);
	x(J) = z;
	x = ProjC(x,r);
    
	stop_criterion1= (norm(A*x-rho_iter*B*x)/(norm_A+abs(rho_iter)*norm_B) < tol1);
	stop_criterion2 = (abs(rho_iter-rho_iter_old) < reltol);    
    stop_criterion = stop_criterion1 || stop_criterion2;
    rho_iter_old = rho_iter;
end

J = (x~=0);
[z,rho_iter]=leadingEig(A(J,J),B(J,J));
z = z/norm(z,2);
x = zeros(Dimension,1);
x(J) = z;
G = 1/rho_iter;

end





function xc = ProjC(x,r)
[~,idx]=maxk(abs(x),r);
xc=sparse(length(x), 1);
xc(idx) = x(idx);
xc = xc./norm(xc,2);
end

function [v,lambda] = leadingEig(A,B)
%% compute the leading eig value and vector of the matrix pair (A,B)
% v is the leading eig vector, lambda is the leading eig value

    %temp = condest(B);
    %if temp> 1e+5 && size(B,1)>5 
    %    disp(temp);
    %end
    [V,D] = eig(A,B,'vector');
    [lambda,temp] = max(abs(D));
    v = V(:,temp);
end

function J = ill_condition_process(B,J,tol3)
	[~,tempR,tempP] = qr(B(J,J),'vector');
	delete_index = tempP(abs(diag(tempR))<tol3*abs(tempR(1,1)));
	if ~isempty(delete_index)
		delete_index2 = find(J,max(delete_index));
		J(delete_index2(delete_index)) = false;
	end
end



function [Q,H] = Arnoldi0(A,k_tot,q1)
% Copyright (c) 2016, Xose Manuel Carreira
% All rights reserved.

n = length(A);
Q = zeros(n,k_tot); 
q1 = q1/norm(q1);
Q(:,1) = q1;
H = zeros(min(k_tot+1,k_tot),n);


for k=1:k_tot
    z = A*Q(:,k);
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        z = z - H(i,k)*Q(:,i);
    end
    if k < n
       H(k+1,k) = norm(z);
       if H(k+1,k) == 0, return, end
       Q(:,k+1) = z/H(k+1,k);
   end
end
end