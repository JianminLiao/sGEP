function [A,B,v] = CCAproblem(d,n,s)
%sparsity s
% n=400; %samples
% d=1000; %total variables num
lambda1=0.9;
% lambda_rest=0.1;

% Blocks of the matrix
d_Block = d/10;
Block = zeros(d_Block);
t = 0.8.^(0:d_Block-1);
for i = 1:d_Block
    Block(i,:) = [t(i:-1:2),t(1:d_Block-i+1)];
end

%real covariance matrix of x and y
Zero = zeros(d_Block);
Sigma_x = [ Block,Zero,Zero,Zero,Zero;
          Zero,Block,Zero,Zero,Zero;
          Zero,Zero,Block,Zero,Zero;
          Zero,Zero,Zero,Block,Zero;
          Zero,Zero,Zero,Zero,Block];
%Sigma xy
vx=ones(d/2,1);
vx(s+1:end)=0;
vx=vx(randperm(d/2));
vx=vx/sqrt(vx'*Sigma_x*vx);

vy=ones(d/2,1);
vy(s+1:end)=0;
vy=vy(randperm(d/2));
vy=vy/sqrt(vy'*Sigma_x*vy);

v=[vx;vy];
v=v/norm(v);

Sigma_xy=lambda1*Sigma_x*vx*vy'*Sigma_x;

Sigma=[Sigma_x,Sigma_xy;Sigma_xy',Sigma_x];
%data matrix
XY=mvnrnd(zeros(d,1),Sigma,n);
X=XY(:,1:d/2);
Y=XY(:,d/2+1:end);
%sample covariance matrix
hat_Sigma=cov([X,Y]);
hat_Sigma_xy=hat_Sigma(1:d/2,d/2+1:end);
hat_Sigma_x=hat_Sigma(1:d/2,1:d/2);
hat_Sigma_y=hat_Sigma(d/2+1:end,d/2+1:end);

A=[zeros(d/2,d/2) hat_Sigma_xy;hat_Sigma_xy' zeros(d/2,d/2)];
B=[hat_Sigma_x zeros(d/2,d/2);zeros(d/2,d/2) hat_Sigma_y];
end