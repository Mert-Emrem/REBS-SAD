% 3 oct
clear all
close all
clc

% LAB 3

%%
% es 2

Ix = 0.0504; 
Iy = 0.0504;
Iz = 0.0109;

I = [Ix 0 0;
     0 Iy 0;
     0 0 Iz];

I_inv = inv(I);

% initial condition
om_x0 = 0.45;
om_y0 = 0.52;
om_z0 = 0.55;

om_0 = [om_x0; om_y0; om_z0];

% analitical sol
lambda= (Iz - Ix)/Ix * om_z0;
om_x_fun = @(t) om_x0* cos(lambda*t) - om_x0* sin(lambda*t);
om_y_fun = @(t) om_y0* sin(lambda*t) - om_y0* cos(lambda*t);
om_z_fun = @(t) om_z0 + 0*t;

% numerical sol
out = sim("LAB_03_3_OCT_SIM_EUL_VECT.slx");
time = out.time;
%
figure(1)
plot(time, om_x_fun(time), 'r');
hold on
plot(time, out.data.omega(:,1), 'b');
title('omega x');


figure(2)
plot(time, om_y_fun(time), 'r');
hold on
plot(time, out.data.omega(:,2), 'b');
title('omega y');


figure(3)
plot(time, om_z_fun(time), 'r');
hold on
plot(time, out.data.omega(:,3), 'b');
title('omega 3');

%% 

%%
% es 3
Ix = 0.0250;
Iy = 0.055;
Iz = 0.07;

om = 2*pi;
om_x0 = om;
om_y0 = 0;
om_z0 = 0;

om_0 = [om_x0; om_y0; om_z0];


