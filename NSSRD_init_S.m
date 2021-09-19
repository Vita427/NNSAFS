function  S=NSSRD_init_S(fea,m)
%% 此函数为非负谱嵌入矩阵S的初始化函数
[N,dim]=size(fea);%400*1024
options.PCARatio=0.75;%0.75
[eigvector, eigvalue] = PCA(fea,options);%1024*19 399*1使用PCA算法对数据进行降维
fea=fea*eigvector;%400*19
label=litekmeans(fea,m,'MaxIter',100,'Replicates',10);%使用聚类算法对数据进行分类
post=zeros(N,m);%400*40
for i=1:N
    post(i,label(i))=1;%Get a pseudo clustering index matrix 
end
S=post';%40*400