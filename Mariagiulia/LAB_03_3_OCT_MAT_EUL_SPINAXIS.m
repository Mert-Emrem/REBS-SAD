% 3 oct 
clear all
close all
clc

% LAB 3

% es 3
Ix = 0.0250;
Iy = 0.055;
Iz = 0.07;

I = diag([Ix, Iy, Iz]);

om = 2*pi;
om_x0 = 0.1;
om_y0 = om;
om_z0 = 0.1;

om_0 = [om_x0; om_y0; om_z0];

out = sim("LAB_03_3_OCT_SIM_EUL_VECT.slx");
time = out.time.data;

figure
plot3(out.omega.data(:,1), out.omega.data(:,2), out.omega.data(:,3));
xlabel('om_x');
ylabel('om_y');
zlabel('om_z');

figure
plot(time, out.omega.data(:,1), 'r');
hold on
plot(time, out.omega.data(:,2), 'b');
plot(time, out.omega.data(:,3), 'y');