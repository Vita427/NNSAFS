function  S=NSSRD_init_S(fea,m)
%% �˺���Ϊ�Ǹ���Ƕ�����S�ĳ�ʼ������
[N,dim]=size(fea);%400*1024
options.PCARatio=0.75;%0.75
[eigvector, eigvalue] = PCA(fea,options);%1024*19 399*1ʹ��PCA�㷨�����ݽ��н�ά
fea=fea*eigvector;%400*19
label=litekmeans(fea,m,'MaxIter',100,'Replicates',10);%ʹ�þ����㷨�����ݽ��з���
post=zeros(N,m);%400*40
for i=1:N
    post(i,label(i))=1;%Get a pseudo clustering index matrix 
end
S=post';%40*400