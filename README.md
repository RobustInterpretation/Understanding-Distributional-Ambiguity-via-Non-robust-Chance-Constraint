# Understanding-Distributional-Ambiguity-via-Non-robust-Chance-Constraint
Scripts for the paper 'Understanding Distributional Ambiguity via Non-robust Chance Constraint'

## Explanation for the codes:

- **CCO_solver.m**: to solve the CCO problem;

- **DRO_solver.m**: to solve the DRO problem using 2nd order reformulation;

- **order4_solver.m**: to solve the DRO problem using 4th order reformulation;

- **RMC_solver_T.m**: to solve the DRO problem using Robust Monte Carlo method;

- **func_rho_N_epsilon_delta_cvx.m**: to find the equivalent ambiguity radius rho given the pair of (epsilon, delta) with 2nd order reformulation solving the DRO problem and center distribution P0 assumed as multivariate normal distribution;

- **func_rho_T_epsilon_delta_cvx.m**: to find the equivalent ambiguity radius rho given the pair of (epsilon, delta) with 2nd order reformulation solving the DRO problem and center distribution P0 assumed as multivariate t distribution;

- **func_rho_N_epsilon_delta_order4.m**: to find the equivalent ambiguity radius rho given the pair of (epsilon, delta) with 4th order reformulation solving the DRO problem and center distribution P0 assumed as multivariate normal distribution;

- **func_rho_T_epsilon_delta_order4.m**: to find the equivalent ambiguity radius rho given the pair of (epsilon, delta) with 4th order reformulation solving the DRO problem and center distribution P0 assumed as multivariate t distribution;

