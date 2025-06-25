addpath("build_lsqlin\")

% optimim problem:

% Optimization function: 
% fo(u) = u^T * u || (minimizing control usage)

% Constraints: 
% - Control Saturation (same for each actuator
bound = 50; % input 

z_des = [3; 0; 3; 0; 0; 0];

time_des = 5.2; %s
s = floor(time_des/ts);

n_u = size(Bd, 2);  % Number of control inputs
n_x = size(Ad, 1);  % Number of states

Rp = Bd;
for i = 2:s
    Rp = [Bd Ad*Rp];
end

u_keep = pinv(Bd)*(z_des - Ad*z_des); % last control input

C_cost = double(eye(n_u * s));  % Ensure C_cost is of type double
d_cost = double(zeros(n_u * s, 1));  % Ensure d_cost is of type double

ub = double(bound * ones(n_u * s, 1));  % Ensure ub is of type double
lb = double(-bound * ones(n_u * s, 1));  % Ensure lb is of type double

Rp = Bd;
for i = 2:s
    Rp = [Bd Ad*Rp];
end

% Correct Aeq and beq
Aeq = double(Ad*Rp);
beq = double(z_des - Ad^(s + 1) * z_init - Bd*u_keep);

% Solve the constrained least squares problem

options = optimoptions('lsqlin', ...
    'TolFun', 1, ...
    'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point');

[u_start, ~, ~, ~] = lsqlin(C_cost, d_cost, [], [], ...
                                    Aeq, beq, lb, ub, [], options);

disp('Found u_start')

[u_opt, ~, ~, ~] = my_lsqlin(C_cost, d_cost, [], [], ...
                                    Aeq, beq, lb, ub, u_start);
disp('Found u_opt')


%% export to simulink
u_opt = [u_opt; u_keep]; % optimal control, constrained
u_opt_reshaped = reshape(u_opt, 2, []);  
time_vector = (0:ts:(size(u_opt_reshaped, 2)-1)*ts)';  
u_opt_mine = [time_vector, u_opt_reshaped'];  
assignin('base', 'u_opt_mine', u_opt_mine); % esporta nel workspace



