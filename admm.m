% X is the matrix: which is dxn
% y is the vector which is nx1
% lambda is the lasso parameter:
function [result, error,stepRes]=admm(X,y,lambda,rho)

[rows,cols]=size(X);
AtAinv=inv(X'*X+rho*eye(cols));
Atb=X'*y;
xk=zeros(cols,1);
yk=zeros(cols,1);
zk=zeros(cols,1);

yk_plus1=xk;
zk_plus1=zk;
xk_plus1=xk;
lrho=lambda/rho;
stepRes=[];

maxiter=1000;
error=[];
epson=1e-6
for i=1:maxiter
   yk=yk_plus1;
   zk=zk_plus1;
   xk=xk_plus1;
   
   xk_plus1=AtAinv*(Atb+rho*zk-yk);
   
   % soft threshold:
   tmp=xk_plus1+yk/rho;
   for j=1:cols
       if tmp(j)>lrho
           zk_plus1(j)=tmp(j)-lrho;
       elseif tmp(j)<-lrho
           zk_plus1(j)=tmp(j)+lrho;
       else
           zk_plus1(j)=0;
       end
   end
   
   yk_plus1=yk+rho*(xk_plus1-zk_plus1);
   tmpError=(xk-xk_plus1)'*(xk-xk_plus1);
   error=[error,tmpError];
   result=xk_plus1;
   stepRes=[stepRes;result'];
   if tmpError<epson
       break;     
   end
     
end
end