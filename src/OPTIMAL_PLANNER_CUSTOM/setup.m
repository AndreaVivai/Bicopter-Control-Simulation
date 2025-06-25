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
       