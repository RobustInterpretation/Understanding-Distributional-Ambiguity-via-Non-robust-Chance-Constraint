function [stock_N_opt_rho,stock_N_opt_x_CCO,stock_N_opt_v_CCO,stock_N_opt_v_CCO_annualized,...
    stock_N_opt_x_DRO,stock_N_opt_v_DRO,stock_N_opt_v_DRO_annualized] = ...
    func_rho_N_epsilon_delta_cvx(epsilon,delta,bound,str_assetclass)
% assume p0 is multivariate normal, estimate its parameters: [m Sigma nv]
tic;
filename = str_assetclass;
data = load(strcat('data_',filename),'-mat');
data = struct2cell(data); 
data = [data{:}]; 

% fit p0 to m-N
m = mean(data); % row vector
Sigma = cov(data);

% CCO solution space
d = length(m);
kN = -norminv(epsilon);
   
[stock_N_opt_x_CCO,stock_N_opt_v_CCO] = CCO_solver(d,m,bound,kN,Sigma,delta);
stock_N_opt_v_CCO_annualized = (stock_N_opt_v_CCO+1).^252-1;

% DRO solution space
phi = 1;
err = 1e-8;
rho_0_ini = 0;
rho_1_ini = 1;

[stock_N_opt_x_DRO_upper,stock_N_opt_v_DRO_upper] = DRO_solver(d,m,rho_0_ini,phi,Sigma,bound);
[stock_N_opt_x_DRO_lower,stock_N_opt_v_DRO_lower] = DRO_solver(d,m,rho_1_ini,phi,Sigma,bound);

if stock_N_opt_v_CCO < stock_N_opt_v_DRO_lower
    disp('Please increase the initial upper bound of rho: rho_1_ini.');
else if stock_N_opt_v_CCO - stock_N_opt_v_DRO_upper > err
        disp('No intersection between CCO and DRO.');
        stock_N_opt_rho = -1;
    else
        stock_N_opt_v_DRO = stock_N_opt_v_DRO_lower;
        rho = rho_1_ini;
        while abs(stock_N_opt_v_DRO - stock_N_opt_v_CCO) > err
            if rho - (rho_0_ini + rho_1_ini)/2 < 1e-14
                break;
            else
                rho = (rho_0_ini + rho_1_ini)/2;
                [stock_N_opt_x_DRO,stock_N_opt_v_DRO] = DRO_solver(d,m,rho,phi,Sigma,bound);
                if stock_N_opt_v_DRO < stock_N_opt_v_CCO
                    rho_1_ini = rho;
                else
                    rho_0_ini = rho;
                end
            end
        end
        stock_N_opt_v_DRO_annualized = (stock_N_opt_v_DRO+1).^252-1;
        stock_N_opt_rho = rho;
    end
end
toc;

