%% SPACECRAFT ATTITUDE DETERMINATION PROJECT
%
% AIM: Simulate the attitude behaviour and determination of a 6U Cubesat in 
% MEO, using Gyroscopes, Sun Sensors and an Horizon Sensor. Simulate the
% attitude control using reaction wheels.
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
mu_E        = 3.99e14;          %[m^3/s^2]      Gravitational constant of Earth


%% Celestial and orbital mechanics
% Parameters needed in the Simulink model to determine the orbit, the
% behaviour of both the spacecraft and the celestial bodies.

e       = 0;                    %[-]        Eccentricity
a       = 16.371e6;              %[m]        Semi-Major Axis
i       = 0.2;                  %[rad]      Inclination
n       = sqrt(G*M_Earth/a^3);  %[rad/s]    Mean angular velocity
T_orb   = 2*pi*sqrt(a^3/mu_E);  %[s]        Orbital Period

% Initial Conditions
w0      = [0;0;n_SC];            %[rad/s]    Initial angular velocity
%A0 =rand(3,3);
A0      = eye(3,3);             %[rad]      Initial attitude matrix
theta0  = 3/2*pi;               %[rad]      Initial angle on obrit

%% Spacecraft specifics
% Parameters that describes the constitution of the spacecraft, and its
% properties related to mass.

% Structure
I                      = 12/23.342 *diag([0.343,0.344,0.286]);             %[kg*m^2]   Inertia matrix
Versor_Surfaces_Matrix  = [  1 , 0 ,-1 , 0 , 0 , 0 , 0 , 0 , 0 , 0  ;
                             0 , 1 , 0 ,-1 , 0 , 0 , 0 , 0 , 0 , 0  ;
                             0 , 0 , 0 , 0 , 1 ,-1 , 1 ,-1 , 1 ,-1  ];      %[-]        Surfaces' direction versors
rhoD_vector             = ones(10,1)*0.1;                                   %[-]        Specular reflectivity coefficent
rhoS_vector             = ones(10,1)*0.5;                                   %[-]        Diffuse reflectivity coefficent
Surfaces_Matrix         = [2e-2;6e-2;2e-2;6e-2;3e-2;3e-2;12e-2;...
                                                    12e-2;12e-2;12e-2];     %[m^2]      Components' surfaces
COM                     = 0.008187;                                         %[m] Shift of the centre of mass on the z axis
position                = 1e-2*[15,0,-COM;0,5,-COM;-15,0,-COM;0,-5,-COM;0,0,10-COM;...
                                    0,0,-10-COM;0,20,10-COM;0,20,10-COM;0,-20,10-COM;0,-20,10-COM]; %[m]        Components' position

%% Perturbations
% Parameters that describes the perturbations that affect the Spacecraft

% SRP Perturbance
P           = F_e/c;                %[N/m^2]    Solar radiation pressure

disp("SRP block implemented");

% Magnetic Perturbance
Theta_m     = deg2rad(11.5);        %[rad]      Magnetic field inclination
j_B         = [0.01;0.05;0.01];     %[A*m^2]    S/C Magetic momentum
g_10        = -29615*10^-9;         %[-]        Gaussian
g_11        = -1728*10^-9;         %[-]        Gaussian
h_11        = 5186*10^-9;           %[-]        Gaussian
H0          = sqrt(h_11^2+g_11^2+g_10^2);

disp("Magnetic Perturbance block implemented");

% Drag Perturbance

C_d = 2.1;                          %[-]        Drag coefficient- range between 1.5and 2.5
A1 = 0.02;
A2 = A1;
A3 = A1;

%% Sensors
% Parameters that describes the behaviour and teh specifics of the sensors

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
HS.FOV                              = pi;                                   %[rad]  Horizon sensor field of view
HS.sample_time                      = 0.1;                                  %[s]    Gyroscope sampling time 

disp("Horizon sensor block implemented");

% Sun Sensor
SunSensor.Misalignment_matrix              = [1,0,0;0,1,0;0,0,1];                  %[rad]  Misallignement Matrix
SunSensor.FOV                              = 2*pi;                                 %[rad]  Horizon sensor field of view
SunSensor.sample_time                      = 0.1;                                  %[s]    Gyroscope sampling time 

disp("Sun sensor block implemented");

%% Actuators
% Parameters that describes the behaviour and specifics of the actuators

% Gyroscopes (3 axis + diagonal config)
RW.A                        =    [1, 0, 0, 1/sqrt(3);
                                  0, 1, 0, 1/sqrt(3);
                                  0, 0, 1, 1/sqrt(3)]; 

RW.A_pseudo_inv             =    [5/6           , -1/6          , -1/6;
                                 -1/6           ,  5/6          , -1/6;
                                 -1/6           , -1/6          ,  5/6;
                                  1/(2*sqrt(3)) , 1/(2*sqrt(3)) , 1/(2*sqrt(3))];

RW.IC                       = [0;0;0;0];                                            %[N*s] Angular momentum initial conditions
RW.Max_torque_saturation    = 0.020;                                                %[N*m] Reaction wheels maximum torque
RW.Max_ang_mom_saturation   = 120/1000;                                            %[N*s] Reaction wheels maximum angu<lar momentum 

disp("Reaction wheels block implemented");

%% Attitude Determination

alpha_1 = 0.5;
alpha_2 =0.5;

%% Attitude controll

% Detumbling
time_detumbling = 1;
n_SC = 2*pi / T_orb;
omega_detumbling = [0; 0; n_SC];
k_d_detumbling = 0.1;

% Slew
time_slew = time_detumbling;
omega_slew = omega_detumbling;
k_d_slew = 0.00001;
k_p_slew = 0.00001;

% Nadir Pointing
time_nadir_pointing = 1;
omega_nadir = omega_detumbling;
k_d_nadir = 0.1;
k_p_nadir = 0.1;

%% PLOTS

% Plot different perturbations acting on the satellite
GG_perturbation         = vecnorm(out.GG_perturbation,2,2);
SRP_perturbation        = vecnorm(out.momentum_SRP,2,2);
Magnetic_perturbation   = vecnorm(out.Magnetic_Perturbation,2,2);

figure
semilogy(out.time,GG_perturbation);
hold on
semilogy(out.time,SRP_perturbation);
hold on
semilogy(out.time,Magnetic_perturbation);
title('Perturbations acting on the S/C');
xlabel('Time [s]');
ylabel('Momentum [N*m]');
ylim([1e-14 1e-4]);
legend('Gravity Gradient Perturbation','Sun Radiation Pressure Perturbation','Magnetic Perturbation');
grid on;


%% PLOTS
% ------------- IMPORTANT!!!! -------------
% RUN EVERYTHING ABOVE HERE + SIMULATION FIRST
% -----------------------------------------

%Plot orbit
% plot3(out.r_N(:,1),out.r_N(:,2),out.r_N(:,3))

% Get r_N, DCM and time from simulink output
r_N = out.r_N;
S_N = out.S_N;
DCM_struct = out.DCM; % Extract logged data structure
DCM_data = DCM_struct.signals.values; % Extract the 3x3xN matrix
time = out.tout(:);

% Satellite model 3D representation
% Main body:

c = 0.008187; % Vertical CoM offset due to solar panels

% 6 vectors to define a rectangular prism
directions = [1, 0, 0; 0, 1, 0; 0, 0, 1; -1, 0, 0; 0, -1, 0; 0, 0, -1];
lens = [0, 0.1, 0.1, 0, 0.1, 0.1]; % x-lengths
widths = [0.1, 0.3405, 0.3405, 0.1, 0.3405, 0.3405]; % y-lengths
heights = [0.2263, 0.2263, 0, 0.2263, 0.2263, 0]; % z-lengths
areas = [widths(1)*heights(1), widths(2)*heights(2), widths(3)*lens(3),...
         widths(4)*heights(4), widths(5)*heights(5), widths(6)*lens(6)]; % surface areas in m^2
positions = [widths(3), 0, 0;
            0, lens(3), 0;
            0, 0, heights(1);
            -widths(3), 0, 0;
            0, -lens(3), 0;
            0, 0, -heights(1)]/2;
length_width_ratios = [0.2263/0.1, 0.2263/0.3405, 0.3405/0.1, 0.2263/0.1, 0.2263/0.3405, 0.3405/0.1];

% Solar panels:

% 2 vectors to define two solar panels
directions = [directions; 0, 0, 1; 0, 0, 1];
lens = [lens 0.3055, 0.3055]; % x-lengths
widths = [widths, 0.345, 0.345]; % y-lengths
areas = [areas, lens(7)*widths(7), lens(8)*widths(8)]; % surface areas in m^2
heights = [heights, 0, 0]; % z-lengths
positions = [positions;
            0, (lens(3)+widths(7))/2, heights(1)/2;
            0, -(lens(3)+widths(7))/2, heights(1)/2];
positions(:,3) = positions(:,3) - c; % all shifted down by vertical CoM offset
length_width_ratios = [length_width_ratios, 0.3055/0.345, 0.3055/0.345]; 

% Makes the satellite visibly large on the plot, 1e^6 is about right
scaleFactor = 2000000;

plotAndAnimateCubesat(positions, directions, areas, length_width_ratios, ...
                      DCM_data, r_N, S_N, time, scaleFactor, 'planet');

function plotAndAnimateCubesat(positions, directions, areas, length_width_ratios, DCM_data, r_N, S_N, time, scaleFactor, mode)
    % Scale satellite properties
    positions = positions * scaleFactor; % Scale positions
    areas = areas * scaleFactor^2; % Scale quads' areas
    
    S_N = 6.378e9*1.2*S_N/norm(S_N);
    r_N_x = r_N(:,1); r_N_y = r_N(:,2); r_N_z = r_N(:,3);
    S_N_x = S_N(:,1); S_N_y = S_N(:,2); S_N_z = S_N(:,3);

    % Number of quads
    numQuads = size(positions, 1);

    % Precompute local vertices for each quad
    patches = gobjects(1, numQuads); % Store patch handles for updating
    localVerticesCell = cell(1, numQuads);

    % Initialize the figure
    hold on;
    grid on;
    view(3);
    axis equal;

    % Initialize trail plot
    trailHandle = scatter3([], [], [], 10, 'r', 'filled');
    trailX = []; trailY = []; trailZ = [];

    % Render Earth if mode is 'planet'
    if strcmp(mode, 'planet')
        R_Earth = 6.378e6;  % Earth's radius in meters

        % Create the sphere for Earth
        [X, Y, Z] = sphere(50); % 50x50 resolution
        X = X * R_Earth; % Scale to Earth's radius
        Y = Y * R_Earth;
        Z = Z * R_Earth;

        % Create the Earth surface with texture
        earth = surf(X, Y, Z, 'EdgeColor', 'none');
        margin = R_Earth * 3; % Adjust margin to fit Earth properly
        xlim([-margin, margin]);
        ylim([-margin, margin]);
        zlim([-margin, margin]);
        set(gca, 'Color', 'k', 'GridColor', [1 1 1], 'XGrid', 'on', 'YGrid', 'on');

        % Load and apply Earth texture
        earthTexture = imread('EarthTexture.jpg');
        earth.CData = flipud(earthTexture); % Flip to align texture correctly
        earth.FaceColor = 'texturemap';
        earth.FaceAlpha = 1; % Fully opaque Earth

        % Initialize the vectors in planet mode
        sunVector = quiver3(0, 0, 0, S_N_x(1), S_N_y(1), S_N_z(1), ...
                            'Color', 'y', 'LineWidth', 2);
        
        scaleBodyAxes = 2e6; % Body axes scale factor in planet mode
        bodyX = quiver3(0, 0, 0, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % X-axis (red)
        bodyY = quiver3(0, 0, 0, 0, 0, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 2); % Y-axis (green)
        bodyZ = quiver3(0, 0, 0, 0, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 2); % Z-axis (blue)

    else
        % Set dynamic axis limits for satellite-only mode
        margin = max(max(abs([r_N_x; r_N_y; r_N_z]))) * 1.5;
        xlim([-margin, margin]);
        ylim([-margin, margin]);
        zlim([-margin, margin]);

        % Initialize the vectors in satellite mode
        scaleBodyAxes_sat = 1000000; % Body axes scale factor in satellite mode
        nadir_scale = 1e-1;
        bodyX = quiver3(0, 0, 0, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5); % X-axis (red)
        bodyY = quiver3(0, 0, 0, 0, 0, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Y-axis (green)
        bodyZ = quiver3(0, 0, 0, 0, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Z-axis (blue)
        nadir_vec = quiver3(0, 0, 0, 0, 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Z-axis (blue)
    end

    % Create and store a quad for each surface
    for i = 1:numQuads
        % Calculate side lengths from scaled area and ratio
        area = areas(i);
        ratio = length_width_ratios(i);
        length = sqrt(area / ratio);
        width = ratio * length;

        % Define local quad vertices (scaled)
        localVertices = [
            -length / 2, -width / 2, 0;
             length / 2, -width / 2, 0;
             length / 2,  width / 2, 0;
            -length / 2,  width / 2, 0
        ]';

        % Define as homogeneous coordinates for easier transformations
        localVertices = [localVertices; ones(1, 4)]; 
        localVerticesCell{i} = localVertices; % Store for later use

        % Plot the initial quad
        patches(i) = patch('Vertices', localVertices(1:3, :)', ...
                           'Faces', [1 2 3 4], ...
                           'FaceColor', 'c', ...
                           'FaceAlpha', 1);
        hold on;
    end

    numFrames = size(DCM_data, 3); % Number of animation frames
    zoomRadius = max(max(abs(positions))) * 2; % Define zoom radius based on satellite's size

    animspeed = 3;
    starting_second = 7850;

    timestep = time(2)-time(1);
    starting_second = starting_second*1/timestep;
   
    for t = starting_second:animspeed:numFrames
        % Get the satellite's position at the current timestep
        satelliteCenter = [r_N_x(t), r_N_y(t), r_N_z(t)];

        % DCM for current frame
        R = DCM_data(:, :, t);

        if strcmp(mode, 'satellite')
            % Update axis limits dynamically to zoom in on the satellite
            xlim([satelliteCenter(1) - zoomRadius, satelliteCenter(1) + zoomRadius]);
            ylim([satelliteCenter(2) - zoomRadius, satelliteCenter(2) + zoomRadius]);
            zlim([satelliteCenter(3) - zoomRadius, satelliteCenter(3) + zoomRadius]);

            % Body axes in the body frame (scaled for visibility)
            bodyX_dir = scaleBodyAxes_sat * R(:, 1);
            bodyY_dir = scaleBodyAxes_sat * R(:, 2);
            bodyZ_dir = scaleBodyAxes_sat * R(:, 3);

            % Update body axes handles
            set(bodyX, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', bodyX_dir(1), 'VData', bodyX_dir(2), 'WData', bodyX_dir(3));
            set(bodyY, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', bodyY_dir(1), 'VData', bodyY_dir(2), 'WData', bodyY_dir(3));
            set(bodyZ, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', bodyZ_dir(1), 'VData', bodyZ_dir(2), 'WData', bodyZ_dir(3));
            set(nadir_vec, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', -r_N_x(t)*nadir_scale, 'VData', -r_N_y(t)*nadir_scale, 'WData', -r_N_z(t)*nadir_scale);
        end


        % Update the Sun vector for the current timestep
        if strcmp(mode, 'planet')
            set(sunVector, 'UData', S_N_x(t), ...
                           'VData', S_N_y(t), ...
                           'WData', S_N_z(t));
        
            % Body axes in the body frame (scaled for visibility)
            bodyX_dir = scaleBodyAxes * R(:, 1);
            bodyY_dir = scaleBodyAxes * R(:, 2);
            bodyZ_dir = scaleBodyAxes * R(:, 3);

            % Update body axes handles
            set(bodyX, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', bodyX_dir(1), 'VData', bodyX_dir(2), 'WData', bodyX_dir(3));
            set(bodyY, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', bodyY_dir(1), 'VData', bodyY_dir(2), 'WData', bodyY_dir(3));
            set(bodyZ, 'XData', satelliteCenter(1), 'YData', satelliteCenter(2), 'ZData', satelliteCenter(3), ...
                       'UData', bodyZ_dir(1), 'VData', bodyZ_dir(2), 'WData', bodyZ_dir(3));

            % Update trail points
            trailX = [trailX, satelliteCenter(1)];
            trailY = [trailY, satelliteCenter(2)];
            trailZ = [trailZ, satelliteCenter(3)];
            set(trailHandle, 'XData', trailX, 'YData', trailY, 'ZData', trailZ);
        end
        
                % Create or update the annotation with the current time
        currentTime = time(t);
        timeString = sprintf('Time: %.2f seconds', currentTime);
        
        if exist('timeAnnotation', 'var') && isvalid(timeAnnotation)
            % Update the annotation text
            set(timeAnnotation, 'String', timeString);
        else
            % Create the annotation text outside the plot
            timeAnnotation = annotation('textbox', [0.1, 0.9, 0.3, 0.1], ... % Position: x, y, width, height
                                        'String', timeString, ...
                                        'FontSize', 12, 'FontWeight', 'bold', ...
                                        'BackgroundColor', 'black', ...
                                        'Color', 'white', ...
                                        'EdgeColor', 'none', ...
                                        'HorizontalAlignment', 'left');
        end

        for i = 1:numQuads
            % Get the DCM for the current timestep 
            R = DCM_data(:, :, t); % Rotation matrix from DCM_data
    
            % Get the spacecraft position in the inertial frame for the current timestep
            spacecraftPosition = [r_N_x(t); r_N_y(t); r_N_z(t)];
    
            % Quad's relative position in the body frame
            relativePosition = positions(i, :)'; % From `positions`
    
            % Rotate the local quad's vertices to align with the direction vector
            direction = directions(i, :)';
            direction = direction / norm(direction); % Normalize direction vector
            z_axis = direction;
    
            % Compute orthogonal basis vectors for the quad's local frame
            if abs(z_axis(3)) < 0.99
                x_axis = cross([0; 0; 1], z_axis);
            else
                x_axis = cross([1; 0; 0], z_axis);
            end
            x_axis = x_axis / norm(x_axis);
            y_axis = cross(z_axis, x_axis);
    
            % Local frame of every quad
            R_quad = [x_axis, y_axis, z_axis];
    
            % Quad redefined using DCM, now it's on body frame of satellite
            R_combined = R * R_quad;
    
            % Transform the quad's relative position to the global frame
            globalPosition = spacecraftPosition + R * relativePosition;

            % Transformation matrix
            T = [R_combined, globalPosition; 0 0 0 1];
    
            % Transform local vertices to global
            globalVertices = T * localVerticesCell{i};
    
            % Update patch vertices
            patches(i).Vertices = globalVertices(1:3, :)';
        end
        
        pause(0.05); 
    end
end


