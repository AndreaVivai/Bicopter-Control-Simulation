function [Ai, bi] = combine_constraints(A_ineq, b_ineq, lb, ub)
    % Initialize inequality constraints
    Ai = [];
    bi = [];

    % Add provided inequality constraints
    if ~isempty(A_ineq)
        Ai = A_ineq;
        bi = b_ineq;
    end

    n = length(lb);

    % Add lower bound constraints
    if ~isempty(lb)
        Ai = [Ai; 
             -eye(n)];
        bi = [bi; 
             -lb];
    end

    % Add upper bound constraints
    if ~isempty(ub)
        Ai = [Ai; 
              eye(n)];
        bi = [bi; 
              ub];
    end
end
