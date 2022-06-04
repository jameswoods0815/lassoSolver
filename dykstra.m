% X is the matrix: which is dxn
% y is the vector which is nx1
% lambda is the lasso parameter:

function [result, error,stepRes]=dykstra(X,y,lambda)
maxIter=1000;
[rows,cols]=size(X);
theta=y;
X1=[X,-X];
error=[];
% each row is vector;
% the residue has 2*cols dimension, each plane has rows number
epson=1e-6;
residual=zeros(2*cols,rows);
stepRes=[];
for i=1:maxIter
   tmperror=[];
   for j=1:2*cols
       
       last_theta=theta;
       theta=halfplaneProject(last_theta-residual(j,:)',X1(:,j),lambda);
       pre_res=residual(j,:)';
       residual(j,:)=(theta-(last_theta-pre_res))';
       
       errorOneStep=(pre_res'-residual(j,:))*(pre_res'-residual(j,:))';
       tmperror=[tmperror,errorOneStep];
       
   end
   errorNow=tmperror*tmperror';
   error=[error,errorNow];
   if errorNow<epson
         break;
       end
   stepRes=[stepRes;theta']; 
end
result=(X'*X)\(X'*(y-theta));

end