function Z=optimal(F,m,nonzero_num)
Z=[];%sparse PC loading matrix
F_init=F;
rho=[];
for i=1:m
    [out1,out2,out3,v]=FullPathData(F,nonzero_num(i));
    Z=[Z,v];%update Z
    F=F-F*(v*v');%deflation
end