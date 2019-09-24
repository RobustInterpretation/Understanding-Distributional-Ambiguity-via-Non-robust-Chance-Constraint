function [Order4_optx,Order4_optval] = order4_solver(d,mu,rho,bound,r)
% input:
% d: the dimension of strategy x; 
% mu: mean of the random vector;
% rho: the ambiguity radius;
% bound: bound for the strategy x;
% r: data source

% output: optimal solution (Order4_optx), optimal value (Order4_optval)

theta = 3; % parameter of alpha divergence
X0 = [ones(d,1)/d;10;0.9]; % column vector: [x;eta1;eta2]
lb = [bound.*ones(d,1);-inf;0]; 
ub = [(1-(d-1)*bound)*ones(d,1);inf;inf]; 
options_new_4 = optimoptions('fmincon','Display','off','FunValCheck','on','MaxFunEvals',1000000,'MaxIter',1000000);   

f = @(X)func_order4(X,theta,rho,r,mu);
[Order4_optx,fval,Order4_exitflag] = fmincon(f,X0,[],[],[ones(1,d),0,0],1,lb,ub,[],options_new_4);
Order4_optval = -fval;

    

