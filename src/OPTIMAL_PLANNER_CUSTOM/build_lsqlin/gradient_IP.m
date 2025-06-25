function gradient = gradient_IP(x, t_IP, Ai, bi, C_func, d_func)
    
    % Gradient of the objective function
    grad_f = C_func' * (C_func * x - d_func);

    % Gradient of the log barrier function
    m = size(Ai, 1); % num. ineq.
    grad_phi = zeros(size(x));

    for i = 1:m
        denom = Ai(i, :) * x - bi(i) ;
        grad_phi = grad_phi - (1 / denom) * Ai(i, :)';
    end

    gradient = t_IP * grad_f + grad_phi;
end
