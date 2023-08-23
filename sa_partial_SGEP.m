%partial support altering for SGEP
function [al_x,AL_X]=sa_partial_SGEP(ini_x,A,B,change_num)
approx_inf=1e2;%a constant
%find the support of in_x;
supp_ini=find(ini_x);
ze_ini=find(ini_x==0);
card_ini=length(supp_ini);%the cardinality of in_x;
if card_ini~=1
    %initialize change records
    in_x_supp_ini=full(ini_x(supp_ini));
    c_ini=in_x_supp_ini'*A(supp_ini,supp_ini)*in_x_supp_ini;
    f_ini=in_x_supp_ini'*B(supp_ini,supp_ini)*in_x_supp_ini;
    q_1_3=c_ini; % a scalar
    q_2_3=f_ini;
    I=ze_ini;
    DIAG_A=diag(A);
    DIAG_B=diag(B);
    Q_1_1=DIAG_A(ze_ini);%AII
    Q_2_1=DIAG_B(ze_ini);%BII
    al_x=full(ini_x);%init value
    AL_X=zeros(length(ini_x),change_num);
    supp_now=supp_ini;
    [~,smallest_id_vec]=mink(abs(al_x(supp_ini)),change_num);
    Q_1_2=A(ze_ini,supp_ini)*al_x(supp_ini); % a length_I by 1 matrix, each row is [Ay]_i.
    Q_2_2=B(ze_ini,supp_ini)*al_x(supp_ini);
    for k=1:change_num
        %set the k-th smallest element of in_x to zero
        %initialize
        supp=supp_now;
        out_id=smallest_id_vec(k);
        j=supp_ini(out_id);
        id1=j==supp;
        supp(id1)=[];
        %update b,c,e,f        
        q_1_3=q_1_3-(A(j,j)*al_x(j)^2+2*al_x(j)*(A(j,supp)*al_x(supp)));
        q_2_3=q_2_3-(B(j,j)*al_x(j)^2+2*al_x(j)*(B(j,supp)*al_x(supp)));
        Q_1_2=Q_1_2-al_x(j)*A(I,j);
        Q_2_2=Q_2_2-al_x(j)*B(I,j);
        al_x(j)=0;
        %
        length_I=length(I);
        Q_1_3=repmat(q_1_3,[length_I,1]);
        Q_2_3=repmat(q_2_3,[length_I,1]);
        D_Q_all_12=Q_1_1.*Q_2_2-Q_1_2.*Q_2_1;
        D_Q_all_13=Q_1_1.*Q_2_3-Q_1_3.*Q_2_1;
        D_Q_all_23=Q_1_2.*Q_2_3-Q_1_3.*Q_2_2;
        DELTA=D_Q_all_13.^2-4*D_Q_all_12.*D_Q_all_23;
        SQRT_DELTA=sqrt(DELTA);
        D_Q_all_12_e_0=D_Q_all_12==0;
        ROOT=(-D_Q_all_13-SQRT_DELTA)./(2.*D_Q_all_12);
        D_Q_all_13_g_0=D_Q_all_13>0;
        D_Q_all_13_l_0=D_Q_all_13<0;
        NE_D_Q_all_23_over_D_Q_all_13=-D_Q_all_23./D_Q_all_13;
        NE_D_Q_all_23_over_D_Q_all_13(isnan(NE_D_Q_all_23_over_D_Q_all_13))=0;
        ENTRY=ROOT+D_Q_all_12_e_0.*(D_Q_all_13_g_0.*approx_inf+D_Q_all_13_l_0.*NE_D_Q_all_23_over_D_Q_all_13);
        VALUE=(Q_1_1.*ENTRY.^2+2.*Q_1_2.*ENTRY+Q_1_3)./(Q_2_1.*ENTRY.^2+2.*Q_2_2.*ENTRY+Q_2_3);
        [~,in_id]=max(VALUE);
        %out_x
        i=I(in_id);
        al_x(i)=ENTRY(in_id);
        AL_X(:,k)=al_x;
        if k~=change_num
            %update supp_now
            supp_now=sort([supp;i]);
            %update Q
            q_1_3=q_1_3+(A(i,i)*ENTRY(in_id)^2+2*ENTRY(in_id)*(A(i,supp)*al_x(supp)));
            q_2_3=q_2_3+(B(i,i)*ENTRY(in_id)^2+2*ENTRY(in_id)*(B(i,supp)*al_x(supp)));
            Q_1_2=Q_1_2+ENTRY(in_id)*A(I,i);
            Q_2_2=Q_2_2+ENTRY(in_id)*B(I,i);
            %delete No.id row
            Q_1_2(in_id,:)=[];
            Q_2_2(in_id,:)=[];
            %update I
            I(in_id)=[];
            %update Q_1_1,Q_2_1
            Q_1_1(in_id)=[];
            Q_2_1(in_id)=[];
        end
    end
else
    [~,in_id]=max(diag(A));
    al_x=zeros(length(A),1);
    al_x(in_id)=1;
    AL_X=al_x;
end
AL_X=sparse(AL_X);
end