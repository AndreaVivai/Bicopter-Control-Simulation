function f = cost_index_IP(x, t_IP, Ai, bi, C_func, d_func)
    
    % Objective function
    residual = C_func * x - d_func;
    f_obj = 0.5 * norm(residual)^2;

    % Log barrier function
    phi = 0;
    for i = 1:size(Ai, 1)
        denom = - Ai(i, :) * x + bi(i);
        phi = phi - log(denom);
    end

    % Total cost
    f = t_IP * f_obj + phi;
end
