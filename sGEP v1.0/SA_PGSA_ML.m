function [X,f_result] = our_method_SGEP_PGSA_M_strict_mono(A,B,eta, options,out_iter_num,x0,epsilon)
%%   Two stage iterative method for SGEP
%     max (x'*A*x)/(x'*B*x)    subject to x'*B*x ~= 0, ||x||_0 <= k
%
% *** inputs ***
% - eta:                    
% - options:                a structure stores user specified parameters which include:
%    -- cardinality_vec:    an m-dimensional vector stores the cadinality for each sparse loading
%    -- optTol:             optimality tolerance
%    -- maxIter:            max number of iteration in each TPower
% - out_iter_num:           iteration limit of our algorithm
% - warm_start:             warm start strategy (0: do not employ warm start; 1: employ warm start)
% - x0:                     initial vector
% - epsilon:                alter the support if (f-f_old)/f_old > epsilon
%
% *** outputs ***
% - X:                      each column of X is a sparse vector
% - f_result:               each column of f_result is the objective value corresponding to the sparse vector

%% Set Parameters
cardinality=options.cardinality;
optTol=options.optTol;
maxIter=options.maxIter;

%% Initialization
X = [];
dim = size(A,1);
% global J %delete
%% Main loop to extract multiple sparse loadings
options_cur.cardinality = cardinality;
options_cur.optTol = optTol;
options_cur.maxIter = maxIter;

%warm start selection
if nargin<6
    [x, G,~] = SGEP_PGSA_M(A,B,cardinality,optTol,eta);
    f=1/G;
else
    [x, G,~] = SGEP_PGSA_M(A,B,cardinality,optTol,eta,x0);
    f=1/G;
end

x_init=x;
a=min(card(x),dim-card(x));
change_num=a;
X = [X,x];
f_result=f;
max_changes=a;
for i=1:out_iter_num-1
    j=1;
    [x2,out_X]=sa_partial_SGEP(x_init,A,B,max_changes);
    while j<=max_changes
        [x, G,~] = SGEP_PGSA_M(A,B,cardinality,optTol,eta,x2);
        f=1/G;
        %Is the value higher?
        if (f-f_result(i))/f_result(i)>epsilon
            %save result
            f_result=[f_result;f];
            x_init=x;
            X = [X,x];
            max_changes=min([card(x),dim-card(x),change_num-1]);%strictly decreasing
            break;
        else
            %decrease the change_num
            change_num=max(max_changes-j,1);
            x2=out_X(:,change_num);
            j=j+1;
        end
    end
    if j>max_changes
        break;
    end
end
x=x_init;%final sparse pc loading
end

function cardinality=card(x)
cardinality=sum(x~=0);
end