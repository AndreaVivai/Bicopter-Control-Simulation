function [x_star, fval, exitflag, iterations] = ...
IP_method(x, Ae, be, epsilon_e, Ai, bi, t_IP, mu_IP, ...
          epsilon_IP, max_iter_IP, epsilon_dm, max_iter_dm, ...
          alpha, beta, C_func, d_func)
   
    % Initialize variables
    m = size(Ai,1); % num of ineuq. constr
    exitflag = 1;
    iterations = 0;

    if norm(Ae*x - be) > epsilon_e || any(Ai*x - bi > 0)
        error('Initial guess x0 is not feasible.');
    end

    for iter_IP = 1:max_iter_IP

        [x, ~, iter_dm] = descent_method_IP(x, Ae, Ai, bi, t_IP, ...
                                            epsilon_dm, max_iter_dm, ...
                                            alpha, beta, C_func, d_func);
        iterations = iterations + iter_dm;

        x_star = x;
        fval = 0.5 * norm(C_func * x - d_func)^2;

        % Check stopping criterion
        if m / t_IP < epsilon_IP 
            disp('Interior Point method solved the optimization problem!')
            exitflag = 0;
            break;

        else
            t_IP = mu_IP * t_IP; % Increase barrier parameter
        end

        disp('')
        disp('---------------------------------------------------------')
        fprintf('iterations: %d\n',iterations);
        fprintf('fval: %d\n',fval);
        fprintf('m/t_IP: %d\n',m/t_IP);

        
    end

    if iter_IP == max_iter_IP
        disp('max iteration for IP reached!')
        exitflag = 0;
    end

end
