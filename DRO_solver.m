function [x,cvx_optval] = DRO_solver(d,m,rho,phi,Sigma,bound)
% input:
% d: the dimension of strategy x; 
% m: mean of the random vector;
% rho: ambiguity radius;
% phi: second derivative of phi function evaluated at 1;
% Sigma: covariance of the random vector;
% bound: bound for the strategy x;

% output: optimal solution (x); optimal value (cvx_optval)

cvx_begin quiet
    variable x(d)
    maximize(x'*m'-sqrt(2*rho/phi)*norm([sqrtm(Sigma)*x]))
    subject to
        sum(x) == 1
        x >= bound
cvx_end