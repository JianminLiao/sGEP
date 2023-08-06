function [X,f_result] = SA_TPM(A, options,out_iter_num)

%%   Enhancing Truncated Power Method for Sparse PCA by Support Altering
%     max x'*A*x    subject to ||x||=1, ||x||_0 <= k
%
% *** inputs ***
% - A:                      p x p covariance matrix
% - options:                a structure stores user specified parameters which include:
%    -- cardinality_vec:    an m-dimensional vector stores the cardinality for each sparse loading
%    -- optTol:             optimality tolerance
%    -- maxIter:            max number of iteration in each TPower
% - out_iter_num:           iteration limit of our algorithm
% 
% *** outputs ***
% - X:                      each column of X is a sparse vector
% - f_result:               each column of f_result is the objective value corresponding to the sparse vector

%% Set Parameters
cardinality_vec=options.cardinality_vec;
optTol=options.OptTol;
maxIter=options.MaxIter;

%% Initialization
X = [];
dim = size(A,1);
%% Main loop to extract multiple sparse loadings
cardinality = cardinality_vec;
options_cur.cardinality = cardinality;
options_cur.optTol = optTol;
options_cur.maxIter = maxIter;

[x, f] = TPower_modified(A, options_cur);
x_init=x;
a=min(card(x),dim-card(x));
change_num=a;
X = [X,x];
f_result=f;
max_changes=a;
for i=1:out_iter_num-1
    j=1;
    [x2,out_X]=sa_partial(x_init,A,max_changes);
    while j<=max_changes
        [x, f] = TPower_modified(A, options_cur,x2);
        %Is the value higher?
        if f>f_result(i)
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
end

function cardinality=card(x)
cardinality=sum(x~=0);
end