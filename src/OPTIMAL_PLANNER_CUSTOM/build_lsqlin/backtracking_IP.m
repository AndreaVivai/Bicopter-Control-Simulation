function t_optimal = backtracking_IP(alpha, beta, x, Deltax, t_IP, Ai, bi, C_func, d_func)
    
    t = 1;

    while true
        x_new = x + t * Deltax;

        % Check if x_new is within feasible region
        % if any(Ai * x_new - bi >= -1e-8)
        %     t = beta * t;
        %     continue;
        % end

        f_new = cost_index_IP(x_new, t_IP, Ai, bi, C_func, d_func);
        f_current = cost_index_IP(x, t_IP, Ai, bi, C_func, d_func);
        G = gradient_IP(x, t_IP, Ai, bi, C_func, d_func);

        % exit condition
        if f_new <= f_current + alpha * t * (G' * Deltax)
            break;
        else
            t = beta * t;
        end

        % Prevent infinite loop
        if t < 1e-12
            error('Line search failed to find a suitable step size.');
        end

    end
    
    t_optimal = t;
end
