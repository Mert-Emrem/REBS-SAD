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
                      DCM_data, r_N, S_N, time, scaleFactor, 'satellite');

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
        margin = R_Earth * 1.75; % Adjust margin to fit Earth properly
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

    animspeed = 100;
    starting_second = 1000;

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
