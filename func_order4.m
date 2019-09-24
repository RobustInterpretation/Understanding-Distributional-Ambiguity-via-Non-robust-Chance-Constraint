function [obj_order4] = func_order4(X,theta,rho,r,mu)

x = X(1:end-2);
eta1 = X(end-1);
eta2 = X(end);

r_portfolio = r*x;
var = moment(r_portfolio,2);
skew = moment(r_portfolio,3);
kurt = moment(r_portfolio,4);

obj_order4 = ((theta-2)*(2*theta-3)/24*(eta1*eta2)^4 + (theta-2)/6*(eta1*eta2)^3 +...
    (1/2+(theta-2)*(2*theta-3)*eta2^2*var/4)*(eta1*eta2)^2 +...
    ((theta-2)*(2*theta-3)*eta2^3*skew/6+(theta-2)*eta2^2*var/2)*(eta1*eta2) + ...
    (theta-2)*(2*theta-3)*eta2^4*kurt/24 + (theta-2)*eta2^3*skew/6 + eta2^2*var/2+rho)/eta2 - mu*x;