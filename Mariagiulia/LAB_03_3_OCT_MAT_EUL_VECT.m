% 3 Oct

% es 1 - Euler equation
% data

clear all
close all
clc

Ix = 0.0700; 
Iy = 0.055;
Iz = 0.025;

I = [Ix 0 0;
     0 Iy 0;
     0 0 Iz];

I_inv = inv(I);

% initial condition
om_x0 = 0.45;
om_y0 = 0.52;
om_z0 = 0.55;

om_0 = [om_x0; om_y0; om_z0];

out = sim("LAB_03_3_OCT_SIM_EUL_VECT.slx");

omx = out.omega(:,1);
omy = out.omega(:,2);
omz = out.omega(:,3);

omx_dot = out.omega_dot(:,1);
omy_dot = out.omega_dot(:,2);
omz_dot = out.omega_dot(:,3);

h = out.hnorm(1);
T = out.E;

den = Ix*Iy*Iz;

Px = ((Iy-Ix).*(h.^2-2*T.*Iz)+(Iz-Ix).*(h.^2-2*T.*Iy))./den;
Py = ((Iz-Iy).*(h.^2-2*T.*Ix)+(Ix-Iy).*(h.^2-2*T.*Iz))./den;
Pz = ((Ix-Iz).*(h.^2-2*T.*Iy)+(Iy-Iz).*(h.^2-2*T.*Ix))./den;
P = [Px, Py, Pz]

Qx = 2*(Iz-Ix)*(Iy-Ix)*Ix/den;
Qy = 2*(Iz-Iy)*(Ix-Iy)*Iy/den;
Qz = 2*(Ix-Iz)*(Iy-Iz)*Iz/den;
Q = [Qx, Qy, Qz]


figure
plot(omx_dot, omx);
hold on
xlabel('\omega_x');
ylabel('\omega dot_x');
xline(0);
yline(0);
grid on
title('\omega_x and \omega dot_x');

figure
plot(omy, omy_dot);
hold on
xlabel('\omega_y');
ylabel('\omega dot_y');
xline(0);
yline(0);
axis on
grid on
title('\omega_y and \omega dot_y');

figure
plot(omz, omz_dot);
hold on
xlabel('\omega_z');
ylabel('\omega dot_z');
xline(0)
yline(0)
grid on
title('\omega_z and \omega dot_z');

T = out.E(1);

Kinetic_energy_ell = @(wx, wy, wz) wx.^2.*Ix./(2.*T) + wy.^2.*Iy./(2.*T) + wz.^2.*Iz./(2.*T) - 1;
Ang_Mom_ell = @(wx, wy, wz) wx.^2.*Ix^2./(h.^2) + wy.^2.*Iy^2./(h.^2) + wz.^2.*Iz^2./(h.^2) - 1;
Polhode_ell = @(wx, wy, wz) wx.^2.*(Ix.*(Ix/h.^2 - 1/(2.*T))) + wy.^2.*(Iy.*(Iy/h.^2 - 1/(2.*T))) + wz.^2.*(Iz.*(Iz/h.^2 - 1/(2.*T)));

figure
[X, Y, Z] = ellipsoid(0, 0, 0, 2*T/Ix, 2*T/Iy, 2*T/Iz, 100);
surf(X, Y, Z);
%axis equal
hold on
[U, V, W] = ellipsoid(0, 0, 0, h^2/Ix^2, h^2/Iy^2, h^2/Iz^2, 100);
surf(U, V, W);
[A, B, C] = ellipsoid(0, 0, 0, 1/(Ix*(Ix/h^2 - 1/(2*T))), 1/(Iy*(Iy/h^2 - 1/(2*T))), 1/(Iz*(Iz/h^2 - 1/(2*T))), 100);
surf(A, B, C)
% [J, K, L] = meshgrid(omx, omy, omz); 
% %%
% M = Polhode_ell(J, K, L); 
% isosurface(J, K, L, M, 0);

%%
figure
hold on
plot3(omy, omz, Kinetic_energy_ell(omx, omy, omz))
axis equal
%[x, y, z] = meshgrid(linspace(2*T(1)/Ix, 2*T(1)/Ix, 50), linspace(-2*T(1)/Iy, 2*T(1)/Iy, 50), linspace(-2*T(1)/Iz, 2*T(1)/Iz, 50));
%v = Kinetic_energy_ell(x, y, z); % Crea il plot dell'ellissoide 
%figure;
%isosurface(x, y, z, v, 0);
%plot3(out.om_dot(:,1), out.om_dot(:,2), out.om_dot(:,3));

% REMBER: looking at om and om-dot is not sufficient: you always need to
% check physical quantitities that conserved and check angular vel in the
% phase space on order and i wx wy wz to see the polhole

% VERIFIY THIS COSED AND PUT A LABEL AS VERIFIED IF EVERYTHING IS OKAY

% Cross product was not available at free so if you don't find something in
% the simulink library is to develope your own block that performs what you
% need 
% ex create your own to compute kinetic energy and save it