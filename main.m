addpath('./dataSets/');
clc;
clear all;
clear memory;
% load('umist.mat');
%   load('AR10P');
% load('Yale_64x64');
% load('Isolet');
% load('PIE10P');
% load('COIL20');
% load('ATT');
% load('TOX-171');
% load('ORL_64x64');
% load('orlraws10P');
rng('default');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fea=X;
gnd=Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nClusts = length(unique(gnd));
c=nClusts;%低维空间的维数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% k-最近邻算法
options = [];%结构体
options.NeighborMode = 'KNN';
options.k = 5;
options.t = 1e+3;
options.WeightMode = 'Heatkernel';
%% 数据图和特征图的权重矩阵Sy和Sp
Sw = constructW(fea',options);%特征图：行向量为特征，列向量为样本
%% 根据权重矩阵得到对角矩阵
Dw = diag(sum(Sw));%将Wp每一列的和作为对角元素
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NITER=30;
k=5;
m=100;%m为选择特征数 ^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化W,H
% parameter = [20,30,40,50,60,70,80,90];
% T_lambda1 = [1e-5 1e-4 1e-3  1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5];
T_lambda1 = [1e0];
T_lambda2 = [1e-3  1e-2 1e-1 1e0 1e1 1e2 1e3];
% T_lambda2 = [1e0];
T_lambda3 = [1e-3  1e-2 1e-1 1e0 1e1 1e2 1e3];
T_lambda4 = [1e0];
result = zeros(98,8);
flag = -1;
tic
% for i_para = 1:8
for i_lambda1 = 1:1
   for i_lambda2 = 1:7
      for i_lambda3 = 1:7
         for i_lambda4 = 1:1
                
%              m = parameter(i_para);
             alpha1 = T_lambda1(i_lambda1);
             alpha2 = T_lambda2(i_lambda2);
             lamda = T_lambda3(i_lambda3);
             beta = T_lambda4(i_lambda4);
             tic
             X=fea';
             [X_new,obj]=test1(X,c,m,alpha1,alpha2,beta,lamda,Sw,Dw,NITER);


for i=1:40
    label=litekmeans(X_new',nClusts,'MaxIter',100,'Replicates',10);
    newres = bestMap(gnd,label);
    AC = length(find(gnd == newres))/length(gnd);
    MIhat=MutualInfo(gnd,label);
    resualt(i,:)=[AC,MIhat];
end
for j=1:2
    a=resualt(:,j);
    ll=length(a);
    temp=[];
    for i=1:ll
        if i<ll-18
            b=sum(a(i:i+19));
            temp=[temp;b];
        end
    end
    [e,f]=max(temp);
    e=e./20;
    MEAN(j,:)=[e,f];
    STD(j,:)=std(resualt(f:f+19,j));
    rr(:,j)=sort(resualt(:,j));
    BEST(j,:)=rr(end,j);
end

             flag = flag+2;
%              result(flag,1)=m;
%              result(flag,2)=alpha1;
%              result(flag,3)=alpha2;
%              result(flag,4)=beta;
%              result(flag,5)=lamda;
%              result(flag:flag+1,6)=STD;
%              result(flag:flag+1,7:8)=MEAN;
             result(flag,1)=alpha1;
             result(flag,2)=alpha2;
             result(flag,3)=beta;
             result(flag,4)=lamda;
             result(flag:flag+1,5)=STD;
             result(flag:flag+1,6:7)=MEAN;

         end
      end
   end
end
% end
toc
