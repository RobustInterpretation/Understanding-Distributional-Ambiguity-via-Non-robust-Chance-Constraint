function [stock_T_opt_rho,stock_T_opt_x_CCO,stock_T_opt_v_CCO,stock_T_opt_v_CCO_annualized,...
    stock_T_opt_x_DRO,stock_T_opt_v_DRO,stock_T_opt_v_DRO_annualized] = ...
    func_rho_T_epsilon_delta_cvx(epsilon,delta,bound,str_assetclass)
% assume p0 is multivariate t, estimate its parameters: [m Sigma nv]
tic;
filename = str_assetclass;
data = load(strcat('data_',filename),'-mat');
data = struct2cell(data); 
data = [data{:}]; %[max(max(data)) min(min(data))]

% fit p0 to m-t
for i=1:size(data,2)
pdt = fitdist(data(:,i),'tlocationscale');
df(i) = pdt.nu;
end
[m Sigma nv] = fitt(data,mean(df)); % m is a row vector

% CCO solution space
d = length(m);
g = (gamma((nv+1)/2)*nv^(nv/2-1))/(gamma(nv/2)*sqrt(pi));
kT = (g./epsilon).^(1/nv);
   
[stock_T_opt_x_CCO,stock_T_opt_v_CCO] = CCO_solver(d,m,bound,kT,Sigma,delta);
stock_T_opt_v_CCO_annualized = (stock_T_opt_v_CCO+1).^252-1;

% DRO solution space
phi = 1;
err = 1e-8;
rho_0_ini = 0;
rho_1_ini = 1;

[stock_T_opt_x_DRO_upper,stock_T_opt_v_DRO_upper] = DRO_solver(d,m,rho_0_ini,phi,Sigma,bound);
[stock_T_opt_x_DRO_lower,stock_T_opt_v_DRO_lower] = DRO_solver(d,m,rho_1_ini,phi,Sigma,bound);

if stock_T_opt_v_CCO < stock_T_opt_v_DRO_lower
    if stock_T_opt_v_CCO == -Inf
        stock_T_opt_rho = -2;
    else
        disp('Please increase the initial upper bound of rho: rho_1_ini.');
    end
else if stock_T_opt_v_CCO > stock_T_opt_v_DRO_upper + err
        disp('No intersection between CCO and DRO.');
        stock_T_opt_rho = -1;
    else
        stock_T_opt_v_DRO = stock_T_opt_v_DRO_lower;
        rho = rho_1_ini;
        while abs(stock_T_opt_v_DRO - stock_T_opt_v_CCO) > err
            if rho - (rho_0_ini + rho_1_ini)/2 < 1e-14
                break;
            else
                rho = (rho_0_ini + rho_1_ini)/2;
                [stock_T_opt_x_DRO,stock_T_opt_v_DRO] = DRO_solver(d,m,rho,phi,Sigma,bound);
                if stock_T_opt_v_DRO < stock_T_opt_v_CCO
                    rho_1_ini = rho;
                else
                    rho_0_ini = rho;
                end
            end
        end
        stock_T_opt_v_DRO_annualized = (stock_T_opt_v_DRO+1).^252-1;
        stock_T_opt_rho = rho;
    end
end
toc;
