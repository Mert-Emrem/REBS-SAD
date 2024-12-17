%% LAB 6
% 7/11/2024

%% es 1

clear all 
close all 
clc

H = 10000;
R = 6.371e+6 + H;
G = 6.67430e-11;
Mt = 5.972e+24;
n2 = G*Mt/R^3;
n = sqrt(n2);

Ix = [0.040; 0.060];
Iy = [0.060; 0.080];
Iz = [0.080; 0.040];

om_x0 = [0; 1e-6];
om_y0 = [0; 1e-6];
om_z0 = [n; n];

tvect = [0: 0.1: 100]';

om_target = [0; 0; n];
A_target_fun = @(t) [cos(n.*t), sin(n.*t), 0.*t; -sin(n.*t), cos(n.*t), 0.*t; 0.*t, 0.*t, 0.*t+1];
A0 = A_target_fun(0);

% case 1
I = diag([Ix(1), Iy(1), Iz(1)]);
om_0 = [om_x0(1); om_y0(1); om_z0(1)];

out = sim('LAB_06_7_NOV_GG_fun.slx');
A_target_m = A_target_fun(out.time);

A_target = [];

for i = 1:length(out.time)
    t = out.time(i);
    A_target_value = A_target_fun(t);
    A_target = cat(3, A_target, A_target_value);
end

error_matrix = [];

plot(out.time, out.error_omega);
for i = 1:length(out.time)
    error_matrix_val = out.A(:,:,i) * A_target(:,:,i)';
    error_matrix = cat(3, error_matrix, error_matrix_val);
end

% FAI IN SUMULINK, NON ASCOLTARE LUDO

%%

clear all 
close all 
clc

H = 10000;
R = 6.371e+6 + H;
G = 6.67430e-11;
Mt = 5.972e+24;
n2 = G*Mt/R^3;
n = sqrt(n2);

Ix = [0.040; 0.060];
Iy = [0.060; 0.080];
Iz = [0.080; 0.040];

om_x0 = [0; 1e-6];
om_y0 = [0; 1e-6];
om_z0 = [n; n];

tvect = [0: 0.1: 100]';
I = diag([Ix(2), Iy(2), Iz(2)]);
om_0 = [om_x0(2); om_y0(2); om_z0(2)];
om_target = [0; 0; n];
A_target_fun = @(t) [cos(n.*t), sin(n.*t), 0.*t; -sin(n.*t), cos(n.*t), 0.*t; 0.*t, 0.*t, 0.*t+1];
A0 = A_target_fun(0);

out = sim('LAB_06_7_NOV_SIM_GGDYN_fun.slx');
out = sim('LAB_06_7_NOV_SIM_GGDYN.slx');





