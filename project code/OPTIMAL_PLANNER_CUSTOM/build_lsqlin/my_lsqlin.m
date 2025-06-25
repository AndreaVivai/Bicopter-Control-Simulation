function [u_opt, resnorm, residual, exitflag, output] = ...
my_lsqlin(C_func, d_func, A_ineq, b_ineq, Aeq, beq, lb, ub, x0)

    % Check if initial guess is provided; else, set to zeros
    if isempty(x0)
        x0 = zeros(size(C_func, 2), 1);
    end

    % Combine inequality constraints with bounds
    [Ai, bi] = combine_constraints(A_ineq, b_ineq, lb, ub);

    % Interior Point Method parameters
    t_IP = 1;          % Initial barrier parameter
    mu_IP = 2;         % Barrier parameter update factor
    epsilon_IP = 1e-9; % Tolerance for the interior point method
    max_iter_IP = 10;  % Maximum iterations for interior point method

    % Descent Method parameters
    alpha = 0.25;       % Backtracking line search parameter
    beta = 0.5;         % Backtracking line search parameter
    epsilon_dm = 1e-9;  % Tolerance for the descent method
    max_iter_dm = 100;  % Maximum iterations for descent method
    t_min = 1e-8;       % Minimum step size

    % Equality constraint tolerance
    epsilon_e = 1e-6;

    [u_opt, ~, exitflag, output] = IP_method(x0, Aeq, beq, epsilon_e, ...
                                             Ai, bi, t_IP, mu_IP, ...
                                             epsilon_IP, max_iter_IP, epsilon_dm, ...
                                             max_iter_dm, alpha, beta, ...
                                             C_func, d_func, t_min);
    % Compute residuals
    residual = C_func * u_opt - d_func;
    resnorm = 0.5 * norm(residual)^2;
end
