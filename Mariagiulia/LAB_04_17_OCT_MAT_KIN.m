%% LAB 4 17 October 


% es 1 - simulate kinematics using the DCM with standard integration
% data
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

A0 = diag([1 1 1]);

om_0 = [om_x0; om_y0; om_z0];

out = sim("LAB_04_17_OCT_SIM_KIN.slx");

%% es_2 quaaternions


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

om_0 = [om_x0; om_y0; om_z0];

out = sim("LAB_04_17_OCT_SIM_KIN_QUAT.slx");


