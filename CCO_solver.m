function [x,cvx_optval] = CCO_solver(d,m,bound,kT,Sigma,delta)
% input:
% d: the dimension of strategy x; 
% m: mean of the random vector;
% bound: bound for the strategy x;
% kT: kappa(epsilon) of t distribution;
% Sigma: covariance of the random vector;
% delta: chance constraint

% output: optimal solution (x); optimal value (cvx_optval)

    cvx_begin quiet
        variable x(d)
        maximize(x'*m')
        subject to
            sum(x) == 1
            x >= bound
            kT*norm([sqrtm(Sigma)*x])-x'*m' <= delta
    cvx_end