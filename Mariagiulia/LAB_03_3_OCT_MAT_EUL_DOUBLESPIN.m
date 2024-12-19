% 3 oct

% LAB 3
clear all
close all
clc
% es 4

Ix = 0.0700;
Iy = 0.055;
Iz = 0.025;
Ir_value = 0.0050;

I = diag([Ix, Iy, Iz]);
Ir = diag([0, 0, Ir_value]);


om_x0 = 1e-6;
om_y0 = 1e-6;
om_z0 = 0.02;
om_r0 = [0, 0, 2*pi];

Tw = [0; 0; 0];

om_0 = [om_x0; om_y0; om_z0];

out = sim("LAB_03_SCALARE.slx");
%time = out.time;

% [row, col] = size(out.omega);
% h = zeros(row,1);
% T = zeros(row,1);

% for i=1:row
%     h(i) = Ix*(out.omega(i,1)) + Iy*(out.omega(i,2) + Iz*(out.omega(i,3)));
%     T(i) = Ix*(out.omega(i,1))^2 + Iy*(out.omega(i,2)^2 + Iz*(out.omega(i,3)))^2;
% end

