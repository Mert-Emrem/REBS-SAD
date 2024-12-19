% 26 sept
% Two orbit dynamics

clear all
close all
clc

%% Data
mu_sun = 1.327e+11;
r0 = [2.029e+8; -1.475e+5; -1.1395e+7];
v0 = [3.021; 24.139; 10.964];

%% Simulation
t0 = 0;
tf = 1000*24*60*60; % days to seconds

out = sim("LAB_02_26_SEPT_SIM_2BODYDYN.slx");

figure
plot3(out.position.data(:,1), out.position.data(:,2), out.position.data(:,3),'r');
title('Orbit');
