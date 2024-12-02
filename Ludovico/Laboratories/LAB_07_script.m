%% SPACECRAFT ATTITUDE DETERMINATION PROJECT
%
% AIM:
%
%
% Made by: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini Mariagiulia
%
% Date:
%
%--------------------------------------------------------------------------

% Clear Memory
clc
clear
close all

%% Constants
% Most important constant needed both in the Simulink model ad to
% calclulate other quantities.

G           = 6.6743e-11;       %[m^3/(kg*s)]   Universal Gravity Constant
R_Earth     = 6.378e6;          %[m]            Earth radius
M_Earth     = 5.972e24;         %[kg]           Earth mass
R_sun       = 149597870000;     %[m]            Sun radius
w_Earth     = 7.29e-5;          %[rad/s]        Earth's angular velocity
T_sun       = 365*24*60*60;     %[s]            Sun Period
n_sun       = 2*pi/T_sun;       %[rad/s]        Sun mean angular velocity
F_e         = 1358;             %[W/m^2]        Solar irradiance
c           = 3e+8;             %[m/s]          Light speed
epsilon     = 0.4093;           %[rad]          Ecliptic inclination [CHECK NAME!!]


%% Celestial and orbital mechanics
% Parameters needed in the Simulink model to determine the orbit, the
% behaviour of both the spacecraft and the celestial bodies.

e       = 0;                    %[-]        Eccentricity
a       = 6.843e6;              %[m]        Semi-Major Axis
i       = 0.1;                  %[rad]      Inclination
n       = sqrt(G*M_Earth/a^3);  %[rad/s]    Mean angular velocity

% Initial Conditions
w0      = [0;0;0.01];           %[rad/s]    Initial angular velocity
A0      = eye(3,3);             %[rad]      Initial attitude matrix
theta0  = 0;                    %[rad]      Initial angle on obrit

%% Spacecraft specifics
% Parameters that describes the constitution of the spacecraft, and its
% properties related to mass.

% Structure
I                       = 10^-2*diag([100.9,25.1,91.6]);                    %[kg*m^2]   Inertia matrix
Versor_Surfaces_Matrix  = [  1 , 0 ,-1 , 0 , 0 , 0 , 1 ,-1 , 1 ,-1  ;
                             0 , 1 , 0 ,-1 , 0 , 0 , 0 , 0 , 0 , 0  ;
                             0 , 0 , 0 , 0 , 1 ,-1 , 0 , 0 , 0 , 0  ];      %[-]        Surfaces' direction versors
rhoD_vector             = ones(10,1)*0.1;                                   %[-]        Specular reflectivity coefficent
rhoS_vector             = ones(10,1)*0.5;                                   %[-]        Diffuse reflectivity coefficent
Surfaces_Matrix         = [6e-2;6e-2;6e-2;6e-2;4e-2;4e-2;12e-2;...
                                                    12e-2;12e-2;12e-2];     %[m^2]      Components' surfaces
position                = 1e-2*[10,0,0;0,10,0;-10,0,0;0,-10,0;0,0,15;...
                                    0,0,-15;0,45,0;0,45,0;0,-45,0;0,-45,0]; %[m]        Components' position

%% Perturbations
% Parameters that describes the perturbations that affect the Spacecraft

% SRP Perturbance
P           = F_e/c;                %[N/m^2]    Solar radiation pressure

disp("SRP block implemented");

% Magnetic Perturbance
Theta_m     = deg2rad(11.5);        %[rad]      Magnetic field inclination
j_B         = [0.01;0.05;0.01];     %[A*m^2]    S/C Magetic momentum
g_01        = -29404.8;             %[-]        Gaussian
g_11        = -1450.9;              %[-]        Gaussian
h_11        = 4652.5;               %[-]        Gaussian

disp("Magnetic Perturbance block implemented");


%% Sensors

% Gyroscope
Gyroscope.A_epsilon                 = [1,0,0;0,1,0;0,0,1];                  %[rad]  Misallignement Matrix
Gyroscope.Gain                      = 1;                                    %[rad]  Gyroscope gain
Gyroscope.sample_time               = 0.1;                                  %[s]    Gyroscope sampling time 
Gyroscope.ARW_standard_deviation    = 0.001;                                %[??]   Gyroscope ARW standard deviation
Gyroscope.ARW_variance              = Gyroscope.ARW_standard_deviation^2;   %[??]   Gyroscope ARW variance
Gyroscope.ARW_mean_value            = 0;                                    %[-]    Gyroscope ARW mean value
Gyroscope.bias_gyro                 = 0;

disp("Gyroscope block implemented");

% Horizon Sensor
HS.Misalignment_matrix              = [1,0,0;0,1,0;0,0,1];                  %[rad]  Misallignement Matrix
HS.FOV                              = 2*pi;                                 %[rad]  Horizon sensor field of view

disp("Horizon sensor block implemented");

%% Attitude Determination
%

alpha_1 = 0.5;
alpha_2 =0.5;

%% PLOTS

%Plot orbit
plot3(out.r_N(:,1),out.r_N(:,2),out.r_N(:,3))