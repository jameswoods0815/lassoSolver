% X is the matrix: which is dxn
% y is the vector which is nx1

function [result,error,stepRes] = coordinateDescent(X, y, lambda)
    X=X';
    [dim,n]=size(X);
    result=zeros(dim,1);
    tmpsum = 2*sum(X.^2, 2);
    tmp = sum(y-X'*result)/n;
    tmpdata = zeros(dim,1);
    epson=1e-8;
    now_error=0;
    iter=0;
    error=[];
    maxIter=1000;
    stepRes=[];
    % coodinate descent
    while(iter<maxIter) 
        last_step_result=result;
        pre_tmp=tmp;
        residual=y-X'*last_step_result-pre_tmp;
        tmp=sum(residual)/n+pre_tmp;
        residual=residual+pre_tmp-tmp;
        
        for k=1:dim
            Xi = X(k,:);
            tmpdata(k,1)=2*Xi*(residual+Xi'*result(k));
            
            % soft threshold
            if(tmpdata(k,1)< -lambda)
                result(k)=(tmpdata(k,1)+lambda)/tmpsum(k);
            elseif(tmpdata(k,1) > lambda)
                result(k)=(tmpdata(k,1)-lambda)/tmpsum(k);
            else
                result(k)=0;
            end
            
            residual = residual - (Xi'*(result(k) - last_step_result(k)));
        end

        last_step_error=now_error;
        now_error =sum((X'*result+tmp-y).^2)+lambda*sum(abs(result));
        error=[error,now_error];
        stepRes=[stepRes;result'];
        if(abs(now_error-last_step_error) <= epson)
            break;
        end     
        iter=iter+1;
        
    end
end