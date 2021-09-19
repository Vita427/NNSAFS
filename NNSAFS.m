function [X_new,obj]=NNSAFS(X,c,m,alpha1,alpha2,beta,lamda,Sw,Dw,NITER)
[d,n]=size(X);
Id=ones(d,1);
Ic=ones(c,1);
W=rand(d,c);
F=orth(rand(n,c));

for iter=1:NITER

%update S
for i=1:n
    for j=1:n
       S(i,j)=exp(-norm(F(i,:)-F(j,:),2)^2/(2*beta));
    end
   S(i,:)=S(i,:)./sum(S(i,:));
end
S=(S+S')./2;
D=diag(sum(S,2));
L=D-S;
%update W
W=W.*((X*F+alpha2*Sw*W)./(X*X'*W+alpha1*Id*Ic'+alpha2*Dw*W+eps));
%update F
% F=F.*((X'*W+miu*F)./(F+2*lamda*L*F+miu*F*F'*F+eps));
A=2*lamda.*L;
B=X'*W;
F=gpi(A,B,1);

tran=0;
for i1=1:n
    for j1=1:n
    tran=tran+S(i1,j1)*log(S(i1,j1));
    end
end

obj(iter)=trace((X'*W-F)*(X'*W-F)')+alpha1*trace(Ic'*W'*Id)+alpha2*trace(W'*(Dw-Sw)*W)+2*lamda*(trace(F'*L*F)+beta*tran);
%     if iter>2
%         if abs(obj(iter)-obj(iter-1))/obj(iter-1)<1e-8
%             break
%         end
%     end
end
score=sum((W.*W),2);
[~,index]=sort(score,'descend');
X_new = X(index(1:m),:);
end