%% LAB 5 24 oct
%

%% es 1

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

om_0 = [om_x0; om_y0; om_z0];
A0 = diag([1 1 1]);
theta0 = pi/4;
phi0 = pi/4;
psi0 = pi/4;

angle0 = [theta0; phi0; psi0];

out = sim("LAB_05_24_OCT_SIM_KIN_COSDIR.slx");
%outA = sim("LAB_04_17_OCT_SIM_KIN.slx");

plot([1:1:1001], out.detA312);
figure
plot([1:1:1001], out.detA313);

flag = zeros(1, 1001);
for i = 1:1:1001
    if out.A312(:,:,i) - out.A313(:,:,i) == 0
       flag(i) = 1;
    end
end


