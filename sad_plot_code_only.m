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
scaleFactor = 750000;

plotAndAnimateCubesat(positions, directions, areas, length_width_ratios, ...
                      DCM_data, r_N, S_N, time, scaleFactor, 'planet');

function plotAndAnimateCubesat(positions, directions, areas, length_width_ratios, DCM_data, r_N, S_N, time, scaleFactor, mode)
    % Scale satellite properties
    positions = positions * scaleFactor; % Scale positions
    areas = areas * scaleFactor^2; % Scale quads' areas
    
    S_N = 6.378e9*0.5*S_N/norm(S_N);
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
        set(gca, 'Color', 'k', 'GridColor', [1 1 1], 'XGrid', 'on', 'YGrid', 'on');

        % Load and apply Earth texture
        earthTexture = imread('EarthTexture.jpg');
        earth.CData = flipud(earthTexture); % Flip to align texture correctly
        earth.FaceColor = 'texturemap';
        earth.FaceAlpha = 1; % Fully opaque Earth

        % Initialize the Sun vector plot (arrow from Earth's center to S_N)
        sunVector = quiver3(0, 0, 0, S_N_x(1), S_N_y(1), S_N_z(1), ...
                            'Color', 'y', 'LineWidth', 2);
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

    numFrames = size(DCM_data, 3); % Number of animation frames
    zoomRadius = max(max(abs(positions))) * 2; % Define zoom radius based on satellite's size

    animspeed = 1;

    for t = 1:animspeed:numFrames
        % Get the satellite's position at the current timestep
        satelliteCenter = [r_N_x(t), r_N_y(t), r_N_z(t)];
        
        if strcmp(mode, 'satellite')
            % Update axis limits dynamically to zoom in on the satellite
            xlim([satelliteCenter(1) - zoomRadius, satelliteCenter(1) + zoomRadius]);
            ylim([satelliteCenter(2) - zoomRadius, satelliteCenter(2) + zoomRadius]);
            zlim([satelliteCenter(3) - zoomRadius, satelliteCenter(3) + zoomRadius]);
        end

        % Update the Sun vector for the current timestep
        if strcmp(mode, 'planet')
            set(sunVector, 'UData', S_N_x(t), ...
                           'VData', S_N_y(t), ...
                           'WData', S_N_z(t));
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
