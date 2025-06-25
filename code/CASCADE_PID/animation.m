% Simulation parameters
scale_factor = 0.01; % Scaling factor for the forces (adjustable)
speed_factor = 0.5;  % Speed factor for the animation (the higher, the faster)

% Extract data 
t = out.t.Time; 
z = out.z.Data;  
u = squeeze(permute(out.u.Data, [3, 1, 2]));  
y_des = out.out_des.Data;  
cost = out.cost.Data;

n_t = length(t);
n_xi = size(z, 1);

% Interpolate u to match the length of t
time_u = linspace(1, n_t, size(u, 1)); % Generate indices over the full range of t
u_interp = interp1(time_u, u, 1:n_t);  % Interpolate u to match the length of t

% Create a static figure with concrete platform
figure;
axis equal;
hold on;

% Set axis limits (invert the y-axis) and remove grid
xlim([-8, 11]);
ylim([-1, 8]);
set(gca, 'Color', [0.8 0.9 1], 'XColor', 'none', 'YColor', 'none'); % No grid or axis lines

% Plot ground
fill([-8 -8 11 11], [-1 0 0 -1], [0 0.5 0], 'EdgeColor', 'none'); % Dark green ground

% Concrete platform (you can adjust the position as needed)
platform_x = -6; % Adjust this value to change the platform's x position
platform_y = 0;  % Adjust this value to change the platform's y position
platform_width = 1;
platform_height = 0.5;
rectangle('Position', [platform_x, platform_y, platform_width, platform_height], 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');

% Add background image and tree (coordinates modifiable manually)
background_image = imread('flat_image.png'); % Load image
tree_x = -4;  % Modifiable
tree_y = 0;  % Modifiable

% Animation loop
for i = 2:(10 * speed_factor):n_t
    cla;
    
    % Plot static ground
    fill([-8 -8 11 11], [-1 0 0 -1], [0 0.5 0], 'EdgeColor', 'none');

    % Plot background image
    imagesc([6, 9], [0, 7], flipud(background_image)); % Ensure y=0 is base, x=6
    
    % Plot concrete platform
    rectangle('Position', [platform_x, platform_y, platform_width, platform_height], 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    
    % Draw tree at desired coordinates
    rectangle('Position', [tree_x, tree_y, 0.6, 2.5], 'FaceColor', [0.3, 0.15, 0.05], 'EdgeColor', 'none'); % Trunk
    rectangle('Position', [tree_x - 0.5, tree_y + 1.5, 1.6, 1.5], 'Curvature', [1, 1], 'FaceColor', [0.1, 0.8, 0.1], 'EdgeColor', 'none'); % Tree top
    
    % Plot the drone and forces
    plotDrone(z(1,1,i), z(3,1,i), z(5,1,i), l, 'k');  
    plotForces(z(1,1,i), z(3,1,i), z(5,1,i), u_interp(i,:), l/2, l/2, scale_factor);

    % Plot desired position
    plot(y_des(i,2), y_des(i,1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    title(['Time: ' num2str(floor(t(i))) ' s']);
    
    % Display velocity and forces
    phi = z(5,1,i);
    velocity_y = z(4,1,i);
    force_right = u_interp(i,1);
    force_left = u_interp(i,2);
    currentcost = cost(i);
    
    %text(-7, 7.5, sprintf('phi: %.2f', phi), 'FontSize', 5);
    %text(-7, 6.5, sprintf('Velocity Y: %.2f', velocity_y), 'FontSize', 10);
    text(-7, 7.5, sprintf('Cost:  %.4e', currentcost), 'FontSize', 10);
    text(-7, 7, sprintf('Force Right: %.2f [N]', force_right), 'FontSize', 10);
    text(-7, 6.5, sprintf('Force Left: %.2f [N]', force_left), 'FontSize', 10);
    
    
    pause(0.05 / speed_factor);
end

function h = plotDrone(x, y, theta, l, color)
    w = 0.2 * l;

    % Drone body
    X = [l/4, -l/4, -l/4, l/4];
    Y = [w/4, w/4, -w/4, -w/4];
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    coords = R * [X; Y];
    X_rot = coords(1, :) + x;
    Y_rot = coords(2, :) + y;
    h = fill(X_rot, Y_rot, color, 'EdgeColor', 'none');

    % Propellers (as ellipses)
    propeller_radius_x = 0.12 * l;
    propeller_radius_y = 0.06 * l;
    propeller_pos = [l/4, -l/4; 0, 0];

    for i = 1:2
        prop_pos_rot = R * propeller_pos(:, i);
        prop_x = prop_pos_rot(1) + x;
        prop_y = prop_pos_rot(2) + y;

        % Plot elliptical propellers
        rectangle('Position', [prop_x - propeller_radius_x, prop_y - propeller_radius_y, 2 * propeller_radius_x, 2 * propeller_radius_y], ...
                  'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceColor', [0.8, 0.8, 0.8, 0.7]);
        
        % Plot lights at the ends of the ellipses
        if i == 1
            light_color = 'r'; % Left (red light)
            rectangle('Position', [prop_x + propeller_radius_x , prop_y - 0.025, 0.05, 0.05], ...
                  'Curvature', [1, 1], 'FaceColor', light_color);
        else
            light_color = 'g'; % Right (green light)
            rectangle('Position', [prop_x - propeller_radius_x - 0.05, prop_y - 0.025, 0.05, 0.05], ...
                  'Curvature', [1, 1], 'FaceColor', light_color);
        end
    end

    % Landing feet
    foot_width = 0.05 * l;
    foot_height = 0.1 * l;
    foot_offset = 0.05 * l;
    rectangle('Position', [x - foot_offset - foot_width/2, y - w/2 - foot_height/2, foot_width, foot_height], 'FaceColor', [0.2 0.2 0.2]);
    rectangle('Position', [x + foot_offset - foot_width/2, y - w/2 - foot_height/2, foot_width, foot_height], 'FaceColor', [0.2 0.2 0.2]);
end

function plotForces(x, y, theta, u, ~, l, scale_factor)
    Fr = scale_factor * u(1);
    Fl = scale_factor * u(2);

    xt = x + (l/2) * cos(theta);
    yt = y + (l/2) * sin(theta);
    xb = x - (l/2) * cos(theta);
    yb = y - (l/2) * sin(theta);

    % Inverti il verso delle forze cambiando i segni nei vettori
    quiver(xt, yt, -Fr * sin(theta), Fr * cos(theta), 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'AutoScaleFactor', 0.5);
    quiver(xb, yb, -Fl * sin(theta), Fl * cos(theta), 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'AutoScaleFactor', 0.5);
end