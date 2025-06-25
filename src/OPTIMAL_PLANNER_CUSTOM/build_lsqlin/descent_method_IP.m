function [x, normGL, iter] = descent_method_IP(x, Ae, Ai, bi, t_IP, ...
                                               epsilon, max_iter, ...
                                               alpha, beta, C_func, d_func)

    for iter = 1:max_iter

        G = gradient_IP(x, t_IP, Ai, bi, C_func, d_func);
        H = hessian_IP(x, t_IP, Ai, bi, C_func);
        [Dx, lambda] = newton(G, H, Ae);

        t = backtracking_IP(alpha, beta, x, Dx, t_IP, Ai, bi, C_func, d_func);

        x = x + t * Dx;
        normGL = norm(G + Ae' * lambda); % Compute norm of Lagrangian gradient

        if normGL < epsilon
            disp('Solved!')
            break;
        end
        disp('')
        fprintf('normGL: %d\n',normGL);
        fprintf('t: %d\n',t);

    end

    if iter == max_iter
        warning('Descent method reached maximum number of iterations.');
    end
end
