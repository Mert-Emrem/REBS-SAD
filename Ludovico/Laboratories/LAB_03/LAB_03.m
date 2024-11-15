%% ESERCIZIO 1
clc
clear

%Data
omega_0 = [0.45 ; 0.52 ; 0.55]; %[rad/s]
I= diag([0.07 , 0.055 , 0.025]); %[kg m^2]
euler_angles_0=[1;1;1];

%% ESERCIZIO 2
clc
clear

%Data
w0  = [0.45 ; 0.52 ; 0.55]; %[rad/s]
I   = diag([0.0504 , 0.0504 , 0.0109]); %[kg m^2]
out = sim("Laboratories\LAB_03\LAB3_Euler_eq.slx");

% Analytical solutions
lam             = ( I(3) - I(1) )/ I(1) * w0(3);
w_analitic_x    = @(t) w0(1)*cos(lam*t) - w0(2)*sin(lam*t);
w_analitic_y    = @(t) w0(1)*sin(lam*t) + w0(2)*cos(lam*t);
w_analitic_z    = @(t) w0(3);



%Evaluate analytical solutions and error

w_x = w_analitic_x(out.tout);
w_y = w_analitic_y(out.tout);
w_z = w_analitic_z(out.tout);
ERR_x = w_x' - out.omega.data(:,1);
ERR_Y = w_y' - out.omega.data(:,2);
ERR_Z = w_z' - out.omega.data(:,3);

%Plot

figure (1)
plot(out.tout, w_x, out.tout , out.omega.data(:,1));
grid on;

figure (2)
plot(out.tout,w_y, out.tout,out.omega.data(:,2));
grid on;

figure(3)
plot(w_z,out.tout, out.omega.data(:,3), out.tout);
grid on;

figure(4)
semilogy(out.tout, ERR_x);
grid on;


%% ESERCIZIO 3
clc
clear

% Data
I       = diag([0.025 , 0.055 , 0.07]); %[kg m^2]
w0    =[2*pi , 1 , 0.1]; %[rad/s]
% w0    =[0.1 , 2*pi , 0.1]; %[rad/s]
%w0    =[0.1 , 0.1 , 2*pi]; %[rad/s]
out     = sim("Laboratories\LAB_03\LAB3_Euler_eq.slx");

% Plot
figure(1)
plot(out.tout,out.omega.Data(:,:))

figure(2)
plot3(out.omega.data(:,1),out.omega.data(:,2),out.omega.data(:,3))
grid on;

%% ESERCIZIO 4
clc
clear

Is      = diag([0.07 , 0.055 , 0.025]);     %[kg m^2]
Ir      = [0 0 0; 0 0 0; 0 0 0.005];         %[kg m^2]
ws_0    = [1e-6 , 1e-6 , 0.02];             %[rad/s]
wr_0    = [ 0 , 0 , 2*pi];                  %[rad/s]
Mr      = [ 0 ; 0 ; 0];




