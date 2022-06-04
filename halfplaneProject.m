% the plane is defined as Ax<=lambda
% A is a column vector; n*1;
% p is a clolum vector; n*1;
% b is the bias
% if a piont is in the plane, return p; if not ,return the projection
% vector;
function result=halfplaneProject(p,A,lambda)
result=p-A*max(A'*p-lambda,0)/(A'*A);
end