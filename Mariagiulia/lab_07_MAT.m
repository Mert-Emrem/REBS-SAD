%% LAB 7
% 14 Nov 24

clear all
close all
clc

Ix = 0.04; 
Iy = 0.06;
Iz = 0.08;

I = [Ix 0 0;
     0 Iy 0;
     0 0 Iz];

% initial condition
om_x0 = 0.45;
om_y0 = 0.52;
om_z0 = 0.55;

A0 = diag([1 1 1]);
%om_0 = [om_x0; om_y0; om_z0];

e = 0.2;
a = 6.843e+6;
i = 0.1;
theta_0 = 0;
n_rot = 0.0012;
eps = deg2rad(23.45);
T = 365*24*3600;
n_sun = 2*pi/T;
r_sun = 149597870000;

H = 10000;
R = 6.371e+6 + H;
G = 6.67430e-11;
Mt = 5.972e+24;
n2 = G*Mt/R^3;
n = sqrt(n2);

om_0 = [0; 0; n];


out = sim('LAB_07_SIM_1')

%%

clear all
close all
clc

I = 1e-2*diag([100.9, 25.1, 91.6]);
A0 = diag([1 1 1]);
H = 10000;
R = 6.843e+6;
G = 6.67430e-11;
Mt = 5.972e+24;
n2 = G*Mt/R^3;
n = sqrt(n2);

e = 0.2;
a = 6.843e+6;
i = 0.1;
theta_0 = 0;
n_rot = 0.0012;
eps = deg2rad(23.45);
T = 365*24*3600;
n_sun = 2*pi/T;
r_sun = 149597870000;

om_0 = [0; 0; n];

Fe = 1358; %w/m^2
c = 299792458; %m/2
P = Fe/c;
A_B1 = 6e-2;
A_B2 = 6e-2;
A_B3 = 6e-2;
A_B4 = 6e-2;
A_B5 = 4e-2;
A_B6 = 4e-2;
A_P1 = 12e-2;
A_P2 = 12e-2;
A_P3 = 12e-2;
A_P4 = 12e-2;
rho_sB = 0.5;
rho_sP = 0.1;
rho_d = 0.1;
N_B1 = [1,0,0]';
N_B2 = [0,1,0]';
N_B3 = [-1,0,0]';
N_B4 = [0,-1,0]';
N_B5 = [0,0,1]';
N_B6 = [0,0,-1]';
N_P1 = [1,0,0]';
N_P2 = [-1,0,0]';
N_P3 = [1,0,0]';
N_P4 = [-1,0,0]';
r_B1 = 1e-2*[10,0,0]';
r_B2 = 1e-2*[0,10,0];
r_B3 = 1e-2*[-10,0,0];
r_B4 = 1e-2*[0,-10,0];
r_B5 = 1e-2*[0,0,15];
r_B6 = 1e-2*[0,0,-15];
r_P1 = 1e-2*[0,45,0];
r_P2 = 1e-2*[0,45,0];
r_P3 = 1e-2*[0,-45,0];
r_P4 = 1e-2*[0,-45,0];

P = 15;

out = sim('LAB_07_SIM_1');

%%
somma = out.M_B1 + out.M_B2 + out.M_B3 + out.M_B4 + out.M_B5 + out.M_B6 + out.M_P1 + out.M_P2 + out.M_P3 + out.M_P4
%%
MAT = [];
for i = 1:1:length(out.time)
    F_B1_test = -P*A_B1*dot(out.S_B(i,:)', N_B1)*((1-rho_sB)*out.S_B(i, :)'+(2*rho_sB*dot(out.S_B(i,:)', N_B1)+2/3*rho_d)*N_B1);
    MAT = [MAT, F_B1_test];
end
MAT = MAT';

%%

clear all
close all
clc

I = 1e-2*diag([100.9, 25.1, 91.6]);
A0 = diag([1 1 1]);

H = 10000;
R = 6.378e+6;
G = 6.67430e-11;
Mt = 5.972e+24;


tvect = [0: 0.1: 100]';


e = 0;
a = 6.843e+6;
i = 0.1;
theta_0 = 0;
n_rot = 0.0012;
eps = deg2rad(23.45);
T = 365*24*3600;
n_sun = 2*pi/T;
r_sun = 149597870000;
n2 = G*Mt/a^3;
n = sqrt(n2);

om_target = [0; 0; n];
A_target_fun = @(t) [cos(n.*t), sin(n.*t), 0.*t; -sin(n.*t), cos(n.*t), 0.*t; 0.*t, 0.*t, 0.*t+1];
%A0 = A_target_fun(0);

om_0 = [0; 0; 0.01];


Fe = 1358; %w/m^2
c = 299792458; %m/2
P = Fe/c;
A_B1 = 6e-2;
A_B2 = 6e-2;
A_B3 = 6e-2;
A_B4 = 6e-2;
A_B5 = 4e-2;
A_B6 = 4e-2;
A_P1 = 12e-2;
A_P2 = 12e-2;
A_P3 = 12e-2;
A_P4 = 12e-2;
rho_sB = 0.5;
rho_sP = 0.1;
rho_d = 0.1;
N_B1 = [1,0,0]';
N_B2 = [0,1,0]';
N_B3 = [-1,0,0]';
N_B4 = [0,-1,0]';
N_B5 = [0,0,1]';
N_B6 = [0,0,-1];
N_P1 = [1,0,0];
N_P2 = [-1,0,0];
N_P3 = [1,0,0];
N_P4 = [-1,0,0];
r_B1 = 1e-2*[10,0,0];
r_B2 = 1e-2*[0,10,0];
r_B3 = 1e-2*[-10,0,0];
r_B4 = 1e-2*[0,-10,0];
r_B5 = 1e-2*[0,0,15];
r_B6 = 1e-2*[0,0,-15];
r_P1 = 1e-2*[0,45,0];
r_P2 = 1e-2*[0,45,0];
r_P3 = 1e-2*[0,-45,0];
r_P4 = 1e-2*[0,-45,0];



g_1_0 = -29404.8;
g_1_1 = -1450.9;
h_1_1 = 4652.5;

j_B = [0.01; 0.05; 0.01];
ang = deg2rad(11.5);
rot_E = 7.29e-5;

P= 15; 
out = sim('LAB_07_SIM_GG_SRP');

%%
figure
semilogy(out.time, vecnorm(out.M_GG, 2, 2), 'b');
hold on
semilogy(out.time, vecnorm(out.M_SRP, 2, 2), 'r');
semilogy(out.time,vecnorm(out.M_Magfiel, 2, 2), 'g');
legend('M_{GG}', 'M_{SRP}', 'M_{MagField}');
grid on
title('Norm of perturbances')

figure
subplot(3,1,1)
semilogy(out.time, out.M_GG(:,1), 'b');
hold on
semilogy(out.time, out.M_GG(:,2), 'r');
semilogy(out.time, out.M_GG(:,3), 'g');
title('M_{GG}');

subplot(3,1,2)
semilogy(out.time, out.M_SRP(:,1), 'b');
hold on
semilogy(out.time, out.M_SRP(:,2), 'r');
semilogy(out.time, out.M_SRP(:,3), 'g');
title('M_{SRP}');

subplot(3,1,3)
semilogy(out.time, out.M_Magfiel(:,1), 'b');
hold on
semilogy(out.time, out.M_Magfiel(:,2), 'r');
semilogy(out.time, out.M_Magfiel(:,3), 'g');
title('M_{MagField}');
