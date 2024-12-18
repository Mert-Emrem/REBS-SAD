%% SPACECRAFT ATTITUDE DETERMINATION PROJECT
%
% TODO:
% - Earth not represented correctly in orientation
%    (0 axial tilt, unrotating)
%
% - Base vectors attached to the inertial frame
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

% Constants
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


% Celestial and orbital mechanics
% Parameters needed in the Simulink model to determine the orbit, the
% behaviour of both the spacecraft and the celestial bodies.

e       = 0;                    %[-]        Eccentricity
a       = 6.843e6;              %[m]        Semi-Major Axis
i       = 0.1;                  %[rad]      Inclination
n       = sqrt(G*M_Earth/a^3);  %[rad/s]    Mean angular velocity

% Initial Conditions
w0      = [1;1;1];           %[rad/s]    Initial angular velocity
A0      = eye(3,3);             %[rad]      Initial attitude matrix
theta0  = 0;                    %[rad]      Initial angle on obrit

% Spacecraft specifics
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

% Perturbations
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

% Drag Perturbance

C_d = 2.1;                          %[-]        Drag coefficient- range between 1.5and 2.5
A1 = 0.02;
A2 = A1;
A3 = A1;
% Sensors

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

% Attitude Determination
%

alpha_1 = 0.5;
alpha_2 =0.5;

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

% 6 vectors to define a rectangular prism
directions = [1, 0, 0; 0, 1, 0; 0, 0, 1; -1, 0, 0; 0, -1, 0; 0, 0, -1];
areas = [0.06, 0.06, 0.04, 0.06, 0.06, 0.04]; % surface areas in m^2
lens = [0, 0.2, 0.2, 0, 0.2, 0.2]; % x-lengths
widths = [0.2, 0, 0.2, 0.2, 0, 0.2]; % y-lengths
heights = [0.3, 0.3, 0, 0.3, 0.3, 0]; % z-lengths
positions = [lens(3), 0, 0;
            0, widths(3), 0;
            0, 0, heights(1);
            -lens(3), 0, 0;
            0, -widths(3), 0;
            0, 0, -heights(1)];
length_width_ratios = [1.5, 1.5, 1, 1.5, 1.5, 1];

% Solar panels:

% (Appended to main body)
% May require trial-error on change

% 2 vectors to define two solar panels
directions = [directions; 0, 0, 1; 0, 0, 1];
areas = [areas, 0.12, 0.12]; % surface areas in m^2
lens = [lens 0.15, 0.15]; % x-lengths
widths = [widths, 0.8, 0.8]; % y-lengths
heights = [heights, 0, 0]; % z-lengths
positions = [positions;
            0, widths(3)+widths(7), 0;
            0, -(widths(3)+widths(7)), 0]/2;
length_width_ratios = [length_width_ratios, 15/80, 15/80]; % .15 cm to .85 cm

% Makes the satellite visibly large on the plot, 1e^6 is about right
scaleFactor = 1000000;

plotAndAnimateCubesat(positions, directions, areas, length_width_ratios, ...
                      DCM_data, r_N, S_N, time, scaleFactor, 'satellite');

function plotAndAnimateCubesat(positions, directions, areas, length_width_ratios, DCM_data, r_N, S_N, time, scaleFactor, mode)
    % Scale satellite properties
    positions = positions * scaleFactor; % Scale positions
    areas = areas * scaleFactor^2; % Scale quads' areas

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
        margin = R_Earth * 1.5; % Adjust margin to fit Earth properly
        xlim([-margin, margin]);
        ylim([-margin, margin]);
        zlim([-margin, margin]);

        % Load and apply Earth texture
        earthTexture = imread('EarthTexture.jpg');
        earth.CData = flipud(earthTexture); % Flip to align texture correctly
        earth.FaceColor = 'texturemap';
        earth.FaceAlpha = 1; % Fully opaque Earth
    else
        % Set dynamic axis limits for satellite-only mode
        margin = max(max(abs([r_N_x; r_N_y; r_N_z]))) * 1.5;
        xlim([-margin, margin]);
        ylim([-margin, margin]);
        zlim([-margin, margin]);
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

    % Animate the quads
    numFrames = size(DCM_data, 3); % Number of animation frames
    zoomRadius = max(max(abs(positions))) * 2; % Define zoom radius based on satellite's size

    for t = 1:numFrames
        % Get the satellite's position at the current timestep
        satelliteCenter = [r_N_x(t), r_N_y(t), r_N_z(t)];
        
        if strcmp(mode, 'satellite')
            % Update axis limits dynamically to zoom in on the satellite
            xlim([satelliteCenter(1) - zoomRadius, satelliteCenter(1) + zoomRadius]);
            ylim([satelliteCenter(2) - zoomRadius, satelliteCenter(2) + zoomRadius]);
            zlim([satelliteCenter(3) - zoomRadius, satelliteCenter(3) + zoomRadius]);
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
    
            % Do the same for the Sun later!!!!!!!!!!!!!!!

            % Transformation matrix
            T = [R_combined, globalPosition; 0 0 0 1];
    
            % Transform local vertices to global
            globalVertices = T * localVerticesCell{i};
    
            % Update patch vertices
            patches(i).Vertices = globalVertices(1:3, :)';
        end
        pause(0.05); % Pause for smooth animation
    end
end
