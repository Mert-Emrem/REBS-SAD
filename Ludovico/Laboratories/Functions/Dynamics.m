function Dynamics(Inertia_matrix,omega_0)
% This function prepare the information for the Dynamics Block
% 
% INPUT:
% I         [3x3]       Inertia matrix
% omega_0   [3x1]       Initial velocity of the satellite
%
% Prepared by: Bernasconi Ludovico
% Date: 20/11/2024
%
% ------------------------------------------------------------------------

%Prepare the variables for the simulink
I = Inertia_matrix;
w0 = omega_0;

disp("Dynamic Block implemented");

end