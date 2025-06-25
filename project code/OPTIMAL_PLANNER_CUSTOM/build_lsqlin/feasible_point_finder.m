function [u_opt, residual, exitflag, output] = ...
feasible_point_finder(A_ineq, b_ineq, Aeq, beq, lb, ub, x0, initialcondition_threshold)

    % Check if initial guess respect the equality constraints
    if isempty(x0)
        disp("Provide an initial guess who respects at least the equality constraints");
        return
    elseif(max(abs(Aeq*x0-beq))>initialcondition_threshold)
        disp("The initial guess doesn't respect the equality constraints");
        return
    end

    % Combine inequality constraints with bounds

    [Ai, bi] = combine_constraints(A_ineq, b_ineq, lb, ub);
    
    % Generate an initial value for the variable 

    s0 = max(Ai*x0-bi);

    % Add a value to s0 to have all the inequality const. inactive

    res = 50;

    s0 = s0 + res;

    % Add a column to the constraints matrixes;
    
    Ai = [Ai -ones(size(Ai,1),1)];
    Aeq = [Aeq zeros(size(Aeq,1),1)];
    
    [Xsub, Idx]=independent_cols(Aeq.');

    Aeq = Aeq(Idx,:);
    beq = beq(Idx);
    
    % New initial condition

    x0 = [x0; s0];

    % Define new cost function:

    C_func = [zeros(size(Ai,2)-1,1);1];
    d_func = 0;

    % Interior Point Method parameters
    t_IP = 1;          % Initial barrier parameter
    mu_IP = 2;        % Barrier parameter update factor
    epsilon_IP = 1e-8; % Tolerance for the interior point method
    max_iter_IP = 50;   % Maximum iterations for interior point method

    % Descent Method parameters
    alpha = 0.25;      % Backtracking line search parameter
    beta = 0.5;        % Backtracking line search parameter
    epsilon_dm = 1e-5; % Tolerance for the descent method
    max_iter_dm = 30;   % Maximum iterations for descent method

    % Equality constraint tolerance
    epsilon_e = 1e-4;

    [u_opt, ~, exitflag, output] = IP_method(x0, Aeq, beq, epsilon_e, ...
                                             Ai, bi, t_IP, mu_IP, ...
                                             epsilon_IP, max_iter_IP, epsilon_dm, ...
                                             max_iter_dm, alpha, beta, ...
                                             C_func, d_func);
    % Compute residuals
    residual = C_func.' * u_opt - d_func;
end
