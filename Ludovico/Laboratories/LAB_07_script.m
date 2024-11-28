clc
clear
close all

e=0.2;
a=6.843e6;
theta0 =0;
rot_rate = 0.0012;
i=0.1;
G=6.6743e-11;
R=6.378e6;
Mt = 5.972e24;
n=sqrt(G*Mt/R^3);
r_sun = 149597870000;
T_sun=365*24*60*60;
n_sun = 2*pi/T_sun;
epsilon=deg2rad(23.45);


% Inertia Matrix
I = 10^-2*diag([100.9,25.1,91.6]);

% Initial Conditions
w0 = [0;0;0.01];
A0 = eye(3,3);

% Ex. Table
Versor_Surfaces_Matrix = [  1   ,0  ,-1 ,0  ,0  ,0  ,1  ,-1 ,1  ,-1;
                            0   ,1  ,0  ,-1 ,0  ,0  ,0  ,0  ,0  ,0;
                            0   ,0  ,0  ,0  ,1  ,-1 ,0  ,0  ,0  ,0];
rhoD_vector = ones(10,1)*0.1;
rhoS_vector = ones(10,1)*0.5;
Surfaces_Matrix = [6e-2;6e-2;6e-2;6e-2;4e-2;4e-2;12e-2;12e-2;12e-2;12e-2];
P=1358/300000000;


position =1e-2*[10,0,0;0,10,0;-10,0,0;0,-10,0;0,0,15;0,0,-15;0,45,0;0,45,0;0,-45,0;0,-45,0];

Theta_m = deg2rad(11.5);
R_Earth = 6378000; %[m]
w_Earth = 7.29e-5;

% Variables
g_01 = -29404.8;
g_11 = -1450.9;
h_11 = 4652.5;

j_B = [0.01;0.05;0.01];

Misalignment_matrix= [1,0,0;0,1,0;0,0,1];
Gain_gyroscope = 1;
Sample_gyroscope =0.1;
Bias_gyroscope = 0;
ARW_standard_deviation =0.1;

Gyroscope = Gyroscope(Misalignment_matrix, Gain_gyroscope, Sample_gyroscope,Bias_gyroscope,ARW_standard_deviation);


%% Horizon Sensor

HS.Misalignment_matrix = [1,0,0;0,1,0;0,0,1];
FOV = 1;

%%
sum(out.momentum_SRP(:,:,1000))



%% Plots
