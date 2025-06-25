function x0 = generate_feasible_xo(Aeq, beq, Ai, bi, lb, ub)

    epsilon_e = 1e-6;  % Tolerance for equality constraints

    % Compute a particular solution x_p to Aeq * x_p = beq
    x_p = Aeq \ beq;

    % Compute the null space of Aeq
    Z = null(Aeq);

    % If the null space is empty, x_p is the only solution to equality constraints
    if isempty(Z)
        x0 = x_p;
        % Check if x_p satisfies inequality constraints and bounds
        if any(Ai * x0 - bi > epsilon_e) || any(x0 < lb - epsilon_e) || any(x0 > ub + epsilon_e)
            error('Cannot find a feasible initial point x0. The problem may be infeasible.');
        end
        return;
    end

    % Set up and solve the LP to adjust x_p within null space to satisfy inequalities
    num_z = size(Z, 2);

    % Objective function: minimize 0 (we just need a feasible point)
    f_lp = zeros(num_z, 1);

    % Inequality constraints
    A_lp = [Ai * Z; -Z; Z];
    b_lp = [bi - Ai * x_p; x_p - lb; ub - x_p];

    % Solve the linear program
    options_lp = optimoptions('linprog', 'Display', 'none');
    [z_opt, ~, exitflag_lp] = linprog(f_lp, A_lp, b_lp, [], [], [], [], options_lp);

    if exitflag_lp ~= 1
        error('Failed to find a feasible initial point x0.');
    end

    % Compute the adjusted feasible x0
    x0 = x_p + Z * z_opt;

    % Verify feasibility
    if norm(Aeq * x0 - beq) > epsilon_e || any(Ai * x0 - bi > epsilon_e) ...
            || any(x0 < lb - epsilon_e) || any(x0 > ub + epsilon_e)
        error('x0 does not satisfy the constraints after adjustment.');
    end
end
