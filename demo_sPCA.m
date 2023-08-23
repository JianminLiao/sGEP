%sPCA on Gaussian data sets
clear
warning off

count=3;
rng("default");
candidate_num=5;
step_size=10;
min_card=10;
max_card=100;
row_num=length(min_card:step_size:max_card);
RESULT=zeros(row_num,candidate_num);
TIME=zeros(row_num,candidate_num);
for j=1:count
    fprintf("%d\n",j)
    result=[];
    time=[];
    len=[];
    %generate the data matrix
    p=300; n=500;                         % p, number of samples; n, number of variables
    D=normrnd(0,1,[p,n]);                       % data matrix
    D=D-repmat((mean(D,1)),p,1);        % Centering of the data
    S=D'*D;
    [~,dim]=size(S);
    [d,ix]=sort(diag(S),'descend');S=S(ix,ix);D=D(:,ix);%sort according to diagonal element for PathSPCA
    for card=min_card:step_size:max_card
        %PathSPCA
        tic
        Z6=optimal(D,1,card);
        PathSPCA_time=toc;
        PathSPCA_value=Z6'*S*Z6;
        
        %Parameter for TPower
        cardinality_vec = card;
        options.cardinality_vec = cardinality_vec;
        options.OptTol = 1e-6;
        options.MaxIter = int32(500);
        options.verbose = int32(0);
        options.initType = int32(2);

        %run TPower
        tic
        [TPower_pcloading,TPower_value] = TPower_SPCA(S, options,false);
        TPower_time=toc;

        % run our backward2_partial_strict_mono without warm-start and no cw-max
        tic
        [U,f_SA_TPM] = SA_TPM(S, options,int32(500));
        SA_TPM_time=toc;
        %our result
        U=full(U);
        SA_TPM_pcloading=U(:,end);
        SA_TPM_value=SA_TPM_pcloading'*S*SA_TPM_pcloading;

        %GCW
        tic
        [gcw_x_out,GCW_value]=cwPCA(S,card,'display',false);
        GCW_time=toc;
        GCW_value=GCW_value(end);

        %PCW
        tic
        [pcw_x_out,PCW_value]=cwPCA(S,card,'type','PCW','display',false);
        PCW_time=toc;
        PCW_value=PCW_value(end);

        result=[result;PathSPCA_value,TPower_value,SA_TPM_value,GCW_value,PCW_value];
        time=[time;PathSPCA_time,TPower_time,SA_TPM_time,GCW_time,PCW_time];
    end
    %calculate the proportion of explained variance
    [largest_eigenvector,largest_eigenvalue]=eigs(S,1,'largestabs','Tolerance',1e-20,'MaxIterations',1000);
    result=result/largest_eigenvalue;
    RESULT=RESULT+result;%Proportion of explained variance
    TIME=TIME+time;
end
RESULT=RESULT/count;
TIME=TIME/count;

%plot PEV
x_vec=min_card:step_size:max_card;
figure
hold on
axis([min_card,max_card,min(min(RESULT)),max(max(RESULT))])
xticks(x_vec);
plot(x_vec,RESULT(:,3),'Marker','^','LineWidth',2)
plot(x_vec,RESULT(:,2),'Marker','+','LineWidth',2)
plot(x_vec,RESULT(:,1),'Marker','o','LineWidth',2)
plot(x_vec,RESULT(:,4),'Marker','square','LineWidth',2)
plot(x_vec,RESULT(:,5),'Marker','x','LineWidth',2)
legend('SA\_TPM','TPM','Path','GCW','PCW','Location','southeast');
xlabel('Sparsity');
ylabel('Proportion of explained variance');
box on
set(gcf,'unit','inches','position',[1,1,5,4.5])
hold off

%plot Computation time
figure
semilogy(x_vec,TIME(:,3),'Marker','^','LineWidth',2)
hold on
semilogy(x_vec,TIME(:,2),'Marker','+','LineWidth',2)
semilogy(x_vec,TIME(:,1),'Marker','o','LineWidth',2)
semilogy(x_vec,TIME(:,4),'Marker','square','LineWidth',2)
semilogy(x_vec,TIME(:,5),'Marker','x','LineWidth',2)
legend('SA\_TPM','TPM','Path','GCW','PCW','Location','northwest');
xlabel('Sparsity');
ylabel('Computation time (Seconds)');
box on
axis([min_card,max_card,min(min(TIME)),max(max(TIME))])
xticks(x_vec);
yticks([0.001,0.002,0.003,0.004,0.005,0.007,0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.3,0.4,0.5,0.7,1,2,3,4,5,7,10,20]);
set(gcf,'unit','centimeters','position',[1,1,14,12])
set(gca,'LineWidth',2);
hold off