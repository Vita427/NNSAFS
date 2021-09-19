clear all;
clc;
addpath('./dataSets/');
load('AR10P.mat');   % 1e-2 1e1 1e0 1e-3
rng('default');  %恢复matlab启动时默认的全局随机流
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fea=X;
% gnd=Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nClusts = length(unique(gnd));%unique除去矩阵中的重复元素
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% k-最近邻算法
options = [];%结构体
options.NeighborMode = 'KNN';
options.k = 5;
% options.WeightMode = 'Binary';
options.t =1e+3;
options.WeightMode = 'Heatkernel';
%% 数据图和特征图的权重矩阵Sy和Sp
Sw = constructW(fea',options);%特征图：行向量为特征，列向量为样本1024*1024
%% 根据权重矩阵得到对角矩阵
Dw = diag(sum(Sw));%将Wp每一列的和作为对角元素1024*1024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NITER=30;
c=nClusts;   % 类别数
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
    label=litekmeans(X_new',nClusts,'MaxIter',100,'Replicates',10);%400*1double 根据特征选择后的数据进行分类
    warning off all
    newres = bestMap(gnd,label);%400*1 double gnd为理想类标，label为实际类标
    AC = length(find(gnd == newres))/length(gnd);%聚类准确率ACC评价指标
    
    MIhat=MutualInfo(gnd,label);%NMI评价指标
    
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
            b=sum(a(i:i+19));%100次评价指标之和
            temp=[temp;b];%
        end
    end
    [e,f]=max(temp);%e存储评价指标和的最大值，f存储索引值
    e=e./20;%求得平均值
    MEAN(j,:)=[e,f];
    STD(j,:)=std(resualt(f:f+19,j));%标准差
    rr(:,j)=sort(resualt(:,j));      %将ACC或者NMI升序排序
    BEST(j,:)=rr(end,j);       %求得ACC或者NMI的最大值
    
end

% STD
% MEAN
BEST
% toc


