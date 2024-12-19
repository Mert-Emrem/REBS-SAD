% 26 SEPT
% spring-mass-dumeper
clear all
close all
clc
%% Data
m = 2;
b = 5;
k = 6;
omega = 0.5;

x0 = 1;
v0 = 0;

%% Simulation
t0 = 0;
tf = 200;
out = sim("LAB_02_26_SEPT_SIM_SMD.slx");

figure
plot(out.time.data, out.x.data, 'r')

figure
plot(out.time.data, out.x_dot.data, 'b');

% se il tempo è nel sistema e non a parte, metti .data per recuperarli dai
% grafici
