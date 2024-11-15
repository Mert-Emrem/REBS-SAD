% ESERCITAZIONE 2

%% ES 1 - Pendulum
clc
clear all

% Data
g           = 9.81;     %[m/s^2]
l           = 5;        %[m]
t0          = 0;        %[s]
tf          = 200;      %[s]
theta0      = 0;        %[rad]
theta_dot0  = pi/4;     %[rad/s]
out         = sim("LAB2_Pendulum.slx","StartTime","t0","StopTime","tf");

%Plot
figure 
plot(out.time.Data,out.Theta.Data);
grid on

figure
plot(out.time.Data,out.Theta_dot.Data);
grid on


%% ES 2 - Mass Spring Dumper

clc
clear all

%Data
m       = 2;    %[kg]
b       = 5;    %[]
k       = 6;    %[]
omega   = 0.5;  %[rad/s]
x0      = 1;    %[m]
v0      = 0;    %[m/s]
t0      = 0;    %[s]
tf      = 200;  %[s]

out         = sim("LAB2_Mass_spring_dumper","StartTime","t0","StopTime","tf");

%Plot
figure 
plot(out.time.Data,out.x.Data);
grid on

figure 
plot(out.time.Data,out.x_dot.Data);
grid on

%% ES 3 - Norm of vector

clc
close all

%Data
v=[0.1817 0.6198 -0.7634];

out = sim("LAB2_Norm_of_vector");

out.norm.Data

%% ES 4 - 2 Body Problem

clc
clear

%Data
mu_sun      = 1.327e11;                     %[km^3/s^2]
r0          = [2.029e8; -1.475e5; -1.1395e7]; %[km]
v0          = [3.021; 24.139; 10.964];        %[km/s]
t0          = 0;                            %[s]
n_giorni    = 100;                          %[giorni]
tf          = 60 * 60 * 24 * n_giorni;

%out         = sim("LAB2_2_Body_problem","StartTime","t0","StopTime","tf");

%Still missing the 3D plot

%% ES 5 - 








