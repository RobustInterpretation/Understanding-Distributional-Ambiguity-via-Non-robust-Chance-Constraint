% robust monte carlo method

% input: r (data source); rho (ambiguity radius)
% output: RMC_optx (optimal solution); RMC_optval (optimal value)

theta = 3;
data = -r; % time horizon * no. of assets
lambda0 = 0.02;

[N,m] = size(data);
N_xsample = 1e2;
x_sample = randfixedsum(m,N_xsample,1,bound,1-(m-1)*bound); % no.of assets * N_xsample
series_return = data*x_sample;
N_lambda = 100;
lambda = linspace(0,lambda0,N_lambda); % row vector
N_sample = N_lambda*N_xsample;

% rho;
N_rho = length(rho);
opt_c = zeros(N_lambda,N_xsample);
feval = zeros(N_lambda,N_xsample);
f_obj_partial = feval;
RMC_optval = zeros(1,N_rho);
RMC_optc = zeros(1,N_rho);
RMC_optlambda = zeros(1,N_rho);
RMC_optx = zeros(m,N_rho);

options = optimoptions('lsqnonlin','Display','off','MaxFunEvals',1000000,'MaxIter',1000000);

tic;
parfor k = 1:N_sample
    [i,j] = ind2sub([N_lambda,N_xsample],k);
    [opt_c(k),feval(k)] = lsqnonlin(@(c)(1-sum((lambda(i)*(theta-1)*series_return(:,j)+c).^(1/(theta-1)))/N),1,...
        -lambda(i)*(theta-1)*min(series_return(:,j)),inf,options);
    f_obj_partial(k) = (theta-1)/(N*theta)*((series_return(:,j))'*(lambda(i)*(theta-1)*series_return(:,j)+opt_c(k)).^(1/(theta-1)))...
        +(opt_c(k)-1)/(lambda(i)*theta*(1-theta)); % +rho/lambda(i);
end
toc;

err = 0.04;
idx = (feval<err);

%
for n = 1:N_rho
    f_obj = f_obj_partial + repmat((rho(n)./lambda)',[1,N_xsample]);
    RMC_optval(n) = -min(f_obj(idx));
    [idx_lambda,idx_x] = find(f_obj == -RMC_optval(n));
    RMC_optc(n) = opt_c(idx_lambda,idx_x);
    RMC_optx(:,n) = x_sample(:,idx_x);
    RMC_optlambda(n) = lambda(idx_lambda);
end
%%