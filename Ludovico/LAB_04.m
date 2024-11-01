%% LABORATORY 04
clc
clear
close all
I = diag([0.07,0.055,0.025]);
w0 = [0.45;0.52;0.55];
A0 =eye(3,3);
q0= [0;0;0;1];




%% Exercise 1
clc
clear
close all

% Inertia Matrix
I = diag([0.07,0.055,0.025]);

% Initial Conditions
w0 = [0.45;0.52;0.55];
A0 = eye(3,3);
q0 = [0;0;0;1];
Euler_angles_0 = [0;0;0];

out = sim("Spacecraft Attitude Dynamics\Laboratories\LAB_04\LAB_04_EX_1.slx");

omega_x=out.omega.Data(:,1);
omega_y=out.omega.Data(:,2);
omega_z=out.omega.Data(:,3);

omega_x_dot=out.omega_dot.Data(:,1);
omega_y_dot=out.omega_dot.Data(:,2);
omega_z_dot=out.omega_dot.Data(:,3);

figure
plot(omega_x,omega_y)
title('Omega x - Omega y');
grid on;

figure
plot(omega_y,omega_z)
title('Omega y - Omega z');
grid on;

figure
plot(omega_z,omega_x)
title('Omega z - Omega x');
grid on;


figure
plot(omega_x_dot,omega_x)
title('Omega x dot - Omega x');
grid on;

figure
plot(omega_y_dot,omega_y)
title('Omega y dot - Omega y');
grid on;

figure
plot(omega_z_dot,omega_z)
title('Omega z dot - Omega z');
grid on;






