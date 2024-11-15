 %% Esercizio 3 - Matlab
clc
clear
clear all

epsilon = -pi/4;
b=[cos(epsilon),sin(epsilon)];

tan_alpha = b(1)/b(2);
alpha_rad=atan(tan_alpha);

alpha_grad=rad2deg(alpha_rad);


fprintf("L'angolo in radianti: %3f \n", alpha_rad)
fprintf("L'angolo in gradi: %3f \n", alpha_grad)

%% ESERCIZIO 1 - Simulink
g=9.81;
l=5;