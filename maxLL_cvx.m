function [x,cvx_status] = maxLL_cvx(f,I_max,fgmean,fgstd,verbose)
% Maximize log-likelihood
% Maximize dot product
% Created by: AA on 08/02/2020
% Last Updated: AA on 08/08/2020

n = size(f,1);

if verbose
    cvx_begin
else
    cvx_begin quiet
end
              variable x(n);
%               maximize( exp(-sum(((f'*x)' - fgmean).^2 ./ (fgstd.^2 + 1),2)) ); % log-likelihood
              maximize( (f'*x)' * fgmean' ); % dot product
%                minimize( norm((f'*x)' - fgmean,2) ) % Mean Square Error
               subject to
                  norm([x;-sum(x)],1) <= 2*I_max;
%                   norm([x;-sum(x)],inf) <= I_max;
    cvx_end