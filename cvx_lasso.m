function result=cvx_lasso(A,b,lambda)
addpath cvx
[m,n]=size(A);
cvx_begin
    cvx_quiet(true);
    variable x(n)
    minimize(0.5*sum_square(A*x-b) + lambda*norm(x,1));
cvx_end
result=x;
end