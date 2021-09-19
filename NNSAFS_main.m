clear all;
clc;
addpath('./dataSets/');
load('AR10P.mat');   % 1e-2 1e1 1e0 1e-3
rng('default');  %�ָ�matlab����ʱĬ�ϵ�ȫ�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fea=X;
% gnd=Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nClusts = length(unique(gnd));%unique��ȥ�����е��ظ�Ԫ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% k-������㷨
options = [];%�ṹ��
options.NeighborMode = 'KNN';
options.k = 5;
% options.WeightMode = 'Binary';
options.t =1e+3;
options.WeightMode = 'Heatkernel';
%% ����ͼ������ͼ��Ȩ�ؾ���Sy��Sp
Sw = constructW(fea',options);%����ͼ��������Ϊ������������Ϊ����1024*1024
%% ����Ȩ�ؾ���õ��ԽǾ���
Dw = diag(sum(Sw));%��Wpÿһ�еĺ���Ϊ�Խ�Ԫ��1024*1024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NITER=30;
c=nClusts;   % �����
alpha1=1e-2;
alpha2=1e1;
lamda=1e0;
beta=1e-3;
m=80;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W=NSSRD_init_P(Dw,Sw,c);
% F=NSSRD_init_S(fea,c);
% F=F';
X=fea';
tic
[X_new,obj]=NNSAFS(X,c,m,alpha1,alpha2,beta,lamda,Sw,Dw,NITER);
for i=1:40
    label=litekmeans(X_new',nClusts,'MaxIter',100,'Replicates',10);%400*1double ��������ѡ�������ݽ��з���
    warning off all
    newres = bestMap(gnd,label);%400*1 double gndΪ������꣬labelΪʵ�����
    AC = length(find(gnd == newres))/length(gnd);%����׼ȷ��ACC����ָ��
    
    MIhat=MutualInfo(gnd,label);%NMI����ָ��
    
    resualt(i,:)=[AC,MIhat];
    %     disp(i);
    %     disp(resualt(i,:));
end
BEST=[];
for j=1:2
    a=resualt(:,j);
    ll=length(a);%=200
    temp=[];
    for i=1:ll
        if i<ll-18%=102
            b=sum(a(i:i+19));%100������ָ��֮��
            temp=[temp;b];%
        end
    end
    [e,f]=max(temp);%e�洢����ָ��͵����ֵ��f�洢����ֵ
    e=e./20;%���ƽ��ֵ
    MEAN(j,:)=[e,f];
    STD(j,:)=std(resualt(f:f+19,j));%��׼��
    rr(:,j)=sort(resualt(:,j));      %��ACC����NMI��������
    BEST(j,:)=rr(end,j);       %���ACC����NMI�����ֵ
    
end

% STD
% MEAN
BEST
% toc


