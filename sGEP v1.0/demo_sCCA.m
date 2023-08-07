clear

count=3;
s=8;%each cannonical vector
rng("default");
candidate_num=4;
step_size=8;
min_card=16;
max_card=10*s;
row_num=length(min_card:step_size:max_card);
RESULT=zeros(row_num,candidate_num);
TIME=zeros(row_num,candidate_num);

for j=1:count
    fprintf("%d\n",j)
    result=[];
    time=[];
    if 1 %TRF CCA
        n=1000;
        p=200;
        [A,B,real_leading_vec] = CCAproblem(n,p,s);
    end
    regularization_par=0.01;
    B_regularized=regularization_par*eye(n)+B;
    norm_A=norm(A);
    norm_B=norm(B_regularized);
    for card=min_card:step_size:max_card
        x0=[ones(card,1);zeros(n-card,1)];
        x0=x0(randperm(n));
        x0=x0/norm(x0);
        if 1
            cardinality = card;
            options.cardinality = cardinality;
            options.optTol = 1e-6;
            options.maxIter = int32(5000);
            options.verbose = int32(0);
            options.initType = int32(2);
            eta=0.99*1/eigs(B_regularized,1);
            epsilon=0;
            if 1                %run TRF
                tic
                [x_TRF,~] =TRF_modified(A+B_regularized,B_regularized,eta,options,x0);
                TRF_time=toc;
                x_TRF=ProjC(x_TRF,card);
                TRF_value=(x_TRF'*A*x_TRF)/(x_TRF'*B_regularized*x_TRF);
            end
            if 1                %run PGSA_ML
                tic
                [x_PGSA_ML, G,iter] = SGEP_PGSA_M(A+B_regularized,B_regularized,card,options.optTol,eta,x0);
                PGSA_ML_time=toc;
                x_PGSA_ML=ProjC(x_PGSA_ML,card);
                f_PGSA_ML=(x_PGSA_ML'*A*x_PGSA_ML)/(x_PGSA_ML'*B_regularized*x_PGSA_ML);
            end
            if 1 % run our_method_SGEP_PGSA_M_strict_mono without warm-start
                tic
                [U7,F_SA_PGSA_ML] = SA_PGSA_ML(A+B_regularized,B_regularized,eta, options,int32(500),x0,epsilon);
                SA_PGSA_ML_time=toc;
                U7(:,end)=ProjC(U7(:,end),card);
                f_SA_PGSA_ML=(U7(:,end)'*A*U7(:,end))/(U7(:,end)'*B_regularized*U7(:,end));
            end
            if 1    %IFTRR
                tic
                v=SGEP_IFTRR(A+B_regularized,B_regularized,card,x0,1e-3,100,norm_A,norm_B);
                IFTRR_time=toc;
                v=ProjC(v,card);
                f_IFTRR=(v'*A*v)/(v'*B_regularized*v);
            end
        end
        result=[result;TRF_value,f_PGSA_ML,f_SA_PGSA_ML,f_IFTRR];
        time=[time;TRF_time,PGSA_ML_time,SA_PGSA_ML_time,IFTRR_time];
    end
    RESULT=RESULT+result;%Proportion of explained variance
    TIME=TIME+time;
end
RESULT=RESULT/count;
TIME=TIME/count;

%plot covariance of the canonical variables
x_vec=min_card:step_size:max_card;
hold on
axis([min_card,max_card,min(min(RESULT)),max(max(RESULT))])
xticks(x_vec);
co=[0,0.447000000000000,0.741000000000000;0.850000000000000,0.325000000000000,0.0980000000000000;0.929000000000000,0.694000000000000,0.125000000000000;0.494000000000000,0.184000000000000,0.556000000000000;0.466000000000000,0.674000000000000,0.188000000000000;0.301000000000000,0.745000000000000,0.933000000000000;0.635000000000000,0.0780000000000000,0.184000000000000];
plot(x_vec,RESULT(:,3),'Marker','^','Color',co(1,:),'LineWidth',2)
plot(x_vec,RESULT(:,2),'Marker','+','Color',co(2,:),'LineWidth',2)
plot(x_vec,RESULT(:,1),'Marker','o','Color',co(3,:),'LineWidth',2)
plot(x_vec,RESULT(:,4),'Marker','square','Color',co(4,:),'LineWidth',2)
legend('SA\_PGSA\_ML','PGSA\_ML','rifle','IFTRR','Location','southeast');
xlabel('Sparsity');
ylabel('Covariance of the canonical variables');
box on
set(gcf,'unit','centimeters','position',[10,10,13,10])

%plot computation time
figure
semilogy(x_vec,TIME(:,3),'Marker','^','Color',co(1,:),'LineWidth',2)
hold on
semilogy(x_vec,TIME(:,2),'Marker','+','Color',co(2,:),'LineWidth',2)
semilogy(x_vec,TIME(:,1),'Marker','o','Color',co(3,:),'LineWidth',2)
semilogy(x_vec,TIME(:,4),'Marker','square','Color',co(4,:),'LineWidth',2)
legend('SA\_PGSA\_ML','PGSA\_ML','rifle','IFTRR','Location','southeast');
xlabel('Sparsity');
ylabel('Computation time (Seconds)');
box on
axis([min_card,max_card,min(min(TIME)),max(max(TIME))])
xticks(x_vec);
yticks([0.001,0.002,0.003,0.004,0.005,0.007,0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.3,0.4,0.5,0.7,1,2,3,4,5,7,10,20,30,40,50,70]);
set(gcf,'unit','centimeters','position',[1,1,14,12])
set(gca,'LineWidth',2);
hold off

function [xc,idx] = ProjC(x,r)
[~,idx]=maxk(x,r,'ComparisonMethod','abs');
xc=sparse(length(x), 1);
x_restrict = x(idx);
xc(idx) = x_restrict / norm(x_restrict);
end