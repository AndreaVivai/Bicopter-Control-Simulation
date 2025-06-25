clear, clc, close all

%parameters
ts = 0.02;

% bicopter constants
m = 5; % [kg]
g = 9.81; % [m/s^2]
l = 0.6; % [m]
Ig = 10*(m*l^2)/3; % [kg*m^2]

% air friction
ni_r = 1;
ni_x = 1;
ni_y = 6;

% change control domain
w = [1,      1;
     l/2  , -l/2];

z_init = [0;0;0;0;0;0];
z_eq = z_init;

u_eq = [m*g; 0];

A1 = [0, 1, 0, 0, 0, 0];
A2 = [0, -ni_x*sin(z_eq(5))/m, 0, 0, -u_eq(1)*cos(z_eq(5))/m - ni_x*z_eq(2)*cos(z_eq(5))/m, 0];
A3 = [0, 0, 0, 1, 0, 0];
A4 = [0, 0, 0, -ni_y*cos(z_eq(5))/m, -u_eq(1)*sin(z_eq(5))/m + ni_y*z_eq(4)*sin(z_eq(5))/m, 0];
A5 = [0, 0, 0, 0, 0, 1];
A6 = [0, 0, 0, 0, 0, -ni_r/Ig];

A = [A1; A2; A3; A4; A5; A6];


B = [0, 0;
     sin(z_eq(5))/m, 0;
     0, 0;
     cos(z_eq(5))/m, 0;
     0, 0
     0,  1/Ig];


C = [0, 0, 1, 0, 0, 0;
     1, 0, 0, 0, 0, 0]; 

D = [0 0;
     0 0];

sys = ss(A,B,C,D);

sysd = c2d(sys, ts, 'zoh');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

% Check if reachability is maintained
R = ctrb(Ad,Bd);

if rank(R) < rank(Ad)
    disp('Reachability lost! The rank of R is:')
    disp(rank(R))
end
       

%% Optimal control Design
addpath("build_lsqlin\")

% optimim problem:

% Optimization function: 
% fo(u) = u^T * u || (minimizing control usage)

% Constraints: 
% - Control Saturation (same for each actuator
bound = 80; % input 

z_des = [1; 0; 1; 0; 0; 0];

time_des = 10; %s
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

% Used just to find a feasable initialization point
% Eventually initial point search will be implemented...
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



