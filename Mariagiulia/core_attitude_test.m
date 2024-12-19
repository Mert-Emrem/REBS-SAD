%% TEST CORE ATTITUDE

clear all
close all
clc

Ix = 0.0700; 
Iy = 0.055;
Iz = 0.025;

I = [Ix 0 0;
     0 Iy 0;
     0 0 Iz];

% initial condition
om_x0 = 0.45;
om_y0 = 0.52;
om_z0 = 0.55;

q0 = [0, 0, 0, 1]';
A0 = diag([1 1 1]);
om_0 = [om_x0; om_y0; om_z0];

out = sim('core_attitude.slx');



% Inertia Matrix
I = diag([0.07,0.055,0.025]);

% Initial Conditions
w0 = [0.45;0.52;0.55];
A0 = eye(3,3);
q0 = [0;0;0;1];
Euler_angles_0 = [0;0;0];

out = sim("LAB_04_EX_1");