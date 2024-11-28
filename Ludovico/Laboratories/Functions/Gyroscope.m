function Gyroscope = Gyroscope(Misalignment_matrix, Gain_gyroscope, Sample_gyroscope,Bias_gyroscope,ARW_standard_deviation)
% This function prepare the information for the Gyroscope Block
% 
% INPUT:
% Misalignment_matrix       [3x3]   Misallignement_Matrix
% Gain_gyroscope            [1]     Gain for the giroscope [adim]
% Sample_gyroscope          [1]     Gyroscope's sampling time [s]
% Bias_gyroscope            [1]     Gyroscope's constant bias error [rad/s]
% ARW_standard_deviation    [1]     Gyroscope's ARW standard deviation [??]
% 
%
%
% Prepared by: Bernasconi Ludovico
% Date: 20/11/2024
%
% ------------------------------------------------------------------------

% Rename variables
Gyroscope.A_epsilon = Misalignment_matrix;
Gyroscope.Gain = Gain_gyroscope;
Gyroscope.sample_time = Sample_gyroscope;

%ARW
Gyroscope.ARW_variance = ARW_standard_deviation^2;
Gyroscope.mean_value = 0;

%RRW
Gyroscope.bias_gyro =Bias_gyroscope;

disp("Gyroscope block implemented");

end