function [x,cvx_status] = wls_dp(A,d,max_I,verbose)
% Maximize log-likelihood
% Maximize dot product
% Created by: AA on 08/07/2020
% Last Updated: AA on 08/08/2020

n = size(A,1);

if verbose
    cvx_begin
else
    cvx_begin quiet
end
              variable x(n);
              minimize( norm((A'*x)' - d',2) ); % mean Square Error
               subject to
                  norm([x;-sum(x)],1) <= 2*max_I;
                  norm([x;-sum(x)],inf) <= max_I;
    cvx_end