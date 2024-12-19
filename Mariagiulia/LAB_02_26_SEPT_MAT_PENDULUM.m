% 26 SEPT 
% MAIN PENDULUM 

clear all
close all
clc

%% DATA
g = 9.81; 
l = 5; 
theta0 = 0;
theta_dot0 = pi/4;

%% Simulation 
t0 = 0;        % start of the sim
tf = 200;      % end of the sim

% after writing this, we switched to simulink

out = sim("LAB_02_26_SEPT_SIM_PENDOLO.slx", "StartTime", "t0", "StopTime", "tf");

%% Plot
figure
plot(out.time, out.theta_dot, 'r');
xlabel('time');
ylabel('theta_dot');

figure
plot(out.time, out.theta, 'b');
xlabel('time');
ylabel('theta');

