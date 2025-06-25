clear, clc, close all

%% parameters

ts = 0.05;

% bicopter constants
m = 5; % [kg]
g = 9.81; % [m/p^2]
l = 0.6; % [m]
Ig = 10*(m*l^2)/3; % [kg*m^2]

% air friction
ni_r = 1;
ni_x = 1;
ni_y = 6;

% change control domain
w = [1,      1;
     l/2  , -l/2];


%% linearized system

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

%% System discretization
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

% Consider a desired state
z_des = [3; 0; 3; ...
        0; 0; 0];

time_des = 5.2; %s
p = floor(time_des/ts);

%% Now i want the system to mantain the reached state, to make it possible, the reached state should

n_u = size(Bd, 2);  % Number of control inputs
n_x = size(Ad, 1);  % Number of states

Rp = Bd;
for i = 2:p
    Rp = [Bd Ad*Rp];
end

% Rp = Bd;
% for i = 2:p
%     Rp = [Ad*Rp Bd];
% end

u_opt = pinv(Rp) * (z_des - Ad^p * z_init); % z_des = Ad^p*z_init + Rp*u_opt

u_keep = pinv(Bd)*(z_des - Ad*z_des); % last control input


% OPTIMIZATION SCRIPT: opt_lsqlin.m

% Objective: Minimize ||u||^2 subject to control saturation, slew rate, 
% and state constraints.

% ----------------------------- INITIALIZATION ----------------------------
% Import parameters

% Consider a desired state
z_des = [3; 0; 3; ...
        0; 0; 0];

time_des = 5.2; %s
p = floor(time_des/ts);

n_u = size(Bd, 2);  % Number of control inputs
n_x = size(Ad, 1);  % Number of states

Rp = Bd;
for i = 2:p
    Rp = [Bd Ad*Rp];
end

u_opt = pinv(Rp) * (z_des - Ad^p * z_init); % z_des = Ad^p*z_init + Rp*u_opt

u_keep = pinv(Bd)*(z_des - Ad*z_des); % last control input


% ----------------------------- CONTROL PARAMETERS -------------------------
% Control Saturation (Bounds for each actuator)
ubounds = [50; 15];    % Upper bounds for F and tau
lbounds = [-50; -15];  % Lower bounds for F and tau

% Control Slew Rate (Maximum rate of change per actuator)
slew_rates = [25; 8];   % Maximum slew rate for F and tau
u_prev = [0; 0];        % Initial control input at t=0

% ----------------------------- STATE CONSTRAINTS --------------------------
% Initialize state bounds (infinity for unconstrained states)
state_lbounds = -inf(n_x, 1);
state_ubounds = inf(n_x, 1);

% Set specific state constraints 
state_lbounds(3) = 0;  % y >= 0
state_ubounds(2) = 5;  % x_dot <= 5
state_ubounds(4) = 5;  % y_dot <= 5
state_lbounds(2) = -5;  % x_dot >= -5
state_lbounds(4) = -5;  % y_dot >= -5
state_ubounds(5) = 0.5;  % theta_dot <= 0.5
state_lbounds(5) = -0.5;  % theta_dot >= -0.5
state_ubounds(6) = 0.1;  % theta_dot <= 0.1
state_lbounds(6) = -0.1;  % theta_dot >= -0.1

% ----------------------------- OPTIMIZATION SETTINGS ----------------------
options = optimoptions('lsqlin', ...
    'MaxIterations', 1e4, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-4, ...
    'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point');

%% ---------------------------- EQUALITY CONSTRAINTS ------------------------
% Compute the reachability matrix Rp inv
Rp = Bd;
for i = 2:p
    Rp = [Ad * Rp Bd];
end

Aeq = Ad*Rp;
beq = z_des - Ad^(p+1) * z_init - Bd * u_keep;
disp('Equality constraints prepared.');

%% -------------------------- CONTROL SATURATION ----------------------------
% Repeat bounds for the entire horizon (p steps)
lb = repmat(lbounds, p, 1);  % Lower bounds
ub = repmat(ubounds, p, 1);  % Upper bounds
disp('Control saturation constraints prepared.');

%% -------------------------- SLEW RATE CONSTRAINTS -------------------------
[A_slew, b_slew] = build_slew_rate_constraints(p, n_u, slew_rates, u_prev, u_keep, ts);
disp('Control slew rate constraints prepared.');

%% -------------------------- STATE CONSTRAINTS -----------------------------
Cv = eye(n_x); % All outputs are directly observable
state_lbounds(isinf(state_lbounds)) = -1e5; % Adjust infinities
state_ubounds(isinf(state_ubounds)) = 1e5;

[Hp, Op] = build_state_constraints(p, n_x, n_u, Ad, Bd, Cv);
state_vlb = repmat(state_lbounds, p, 1);
state_vub = repmat(state_ubounds, p, 1);

A_state = [Hp; -Hp];
b_state = [state_vub - Op * z_init; -state_vlb + Op * z_init];
disp('State constraints prepared.');

%% -------------------------- INEQUALITY CONSTRAINTS ------------------------
% Combine all inequality constraints
A_ineq = [A_slew; A_state];
b_ineq = [b_slew; b_state];

% -------------------------- OBJECTIVE FUNCTION ----------------------------
C_func = eye(n_u * p);          % Minimize ||u||^2 (control effort)
d_func = zeros(n_u * p, 1);     

%% ---------------------------- OPTIMIZATION CALL ---------------------------
disp('Optimization start');
[u_opt, resnorm, residual, exitflag, output] = lsqlin(C_func, d_func, ...
                                                      A_ineq, b_ineq, ...
                                                      Aeq, beq, lb, ub, [], options);

if exitflag <= 0
    warning('Optimization did not converge: %s', output.message);
end

%% -------------------------- EXPORT TO SIMULINK ----------------------------
u_reshaped = reshape(u_opt, n_u, p)';
u_reshaped = [u_reshaped; u_keep'];
time_vector = (0:ts:(p) * ts)';
u_opt_c = [time_vector, u_reshaped];
assignin('base', 'u_opt_vanilla', u_opt_c);

%% ------------------------ FUNCTION DEFINITIONS ----------------------------

function [A_slew, b_slew] = build_slew_rate_constraints(p, n_u, slew_rates, u_prev, u_keep, ts)
    % BUILD_SLEW_RATE_CONSTRAINTS: Constructs the slew rate constraints for optimization
    % u_prev    - previous control input
    % u_keep    - control input to maintain the state
    
    A_slew = [];
    b_slew = [];

    for k = 1:p+1
        for i = 1:n_u
            idx = (k - 1) * n_u + i;
            if k == 1
                % Slew rate between u(1) and u_prev
                % positive slew (acceleration)
                A_slew = [A_slew; build_slew_row(idx, 1, n_u * p)];
                b_slew = [b_slew; slew_rates(i) * ts + u_prev(i)];
                % negative slew (deceleration)
                A_slew = [A_slew; build_slew_row(idx, -1, n_u * p)];
                b_slew = [b_slew; slew_rates(i) * ts - u_prev(i)];
            elseif k <= p
                % Slew rate between u(k) and u(k-1)
                A_slew = [A_slew; build_slew_row(idx, 1, n_u * p, idx - n_u, -1)];
                b_slew = [b_slew; slew_rates(i) * ts];
                A_slew = [A_slew; build_slew_row(idx, -1, n_u * p, idx - n_u, 1)];
                b_slew = [b_slew; slew_rates(i) * ts];
            else
                % Slew rate between u(p+1) and u(p)
                A_slew = [A_slew; build_slew_row(idx - n_u, 1, n_u * p)];
                b_slew = [b_slew; slew_rates(i) * ts + u_keep(i)];
                A_slew = [A_slew; build_slew_row(idx - n_u, -1, n_u * p)];
                b_slew = [b_slew; slew_rates(i) * ts - u_keep(i)];
            end
        end
    end
end

function row = build_slew_row(index, sign, row_length, prev_idx, prev_sign)
    % BUILD_SLEW_ROW: Helper function to create a row for A_slew matrix
    % index     - current index for the control
    % sign      - sign of the current control (1 or -1)
    % row_length- total length of the row (n_u * p)
    % prev_idx  - (optional) previous control index for slew rate calculation
    % prev_sign - (optional) sign for the previous control input

    row = zeros(1, row_length);
    row(index) = sign;
    if nargin > 3
        row(prev_idx) = prev_sign;
    end
end

function [Hp, Op] = build_state_constraints(p, n_x, n_u, Ad, Bd, Cv)
    % BUILD_STATE_CONSTRAINTS: Constructs Hp and Op matrices for state constraints

    Hp = zeros(n_x * p, n_u * p); % Initialize Hp matrix
    Op = zeros(n_x * p, n_x);     % Initialize Op matrix

    for i = 1:p
        for j = 1:i
            Hp((i-1)*n_x+1 : i*n_x, (j-1)*n_u+1 : j*n_u) = Cv * Ad^(i-j) * Bd;
        end
        Op((i-1)*n_x+1 : i*n_x, :) = Cv * Ad^i;
    end
end

