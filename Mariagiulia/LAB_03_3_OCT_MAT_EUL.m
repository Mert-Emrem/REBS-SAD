% 3 Oct

% es 1 - Euler equation
% data
Ix = 0.0700; 
Iy = 0.055;
Iz = 0.025;

% initial condition
om_x0 = 0.45;
om_y0 = 0.52;
om_z0 = 0.55;

om_0 = [om_x0, om_y0, om_z0];

out = sim("LAB_03_3_OCT_SIM_EUL.slx");


