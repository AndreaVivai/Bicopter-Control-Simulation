function H = hessian_IP(x, t_IP, Ai, bi, C_func)

    % Hessian of the objective function
    H_f = C_func' * C_func;

    % Hessian of the log barrier function
    m = size(Ai, 1); % num ineq.
    H_phi = zeros(length(x));

    for i = 1:m
        denom = bi(i) - Ai(i, :) * x;
        H_phi = H_phi + (Ai(i, :)' * Ai(i, :)) / denom^2;
    end

    H = t_IP * H_f + H_phi;
end
