%% Wind Turbine Project
clc
clear   

%% Import all Data and Variables
% Turbine Operational Values
power_rated = 2.5e6; % [W] rated turbine power
R = 48; % [m] blade radius
hub_height = 80.4; % [m] hub height
speed_cut_in = 4; % [m/s] cut in wind speed
speed_rated = 11; % [m/s] rated wind speed
speed_cut_out = 25; % [m/s] cut out wind speed
rotor_speed_range = [9.6, 15.5]; % [RPM] min and max rotor speed
n_blades = 3; % number of blades
alpha = 1/3; % axial induction factor [-] 

% Fluid and Material Properties
dens_air = 1.1; % [kg/m^3] air density
visc_air = 1.8e-5; % [N s/m^2] dynamic viscosity of air
dens_steel = 7850; % [kg/m^3] steel density
strength_tens_steel = 450e6; % [Pa] tensile strength of steel
strength_yield_steel = 345e6; % [Pa] yield strength of steel

% Import Turbine Data
Blade_data = readtable('BladeProfile.xlsx'); % Blade data
tower_data = readtable('towerSpecs.xlsx'); % Tower data

% Airfoil Data stored in structure 
airfoils = struct();
airfoils.DU91 = readtable('DU91-W2-250.xlsx');
airfoils.DU93 = readtable('DU93-W-210.xlsx');
airfoils.DU96 = readtable('DU96-W-180.xlsx');
airfoils.DU97 = readtable('DU97-W-300.xlsx');

% Convert mm to m and store blade data 
radius = Blade_data.DistanceFromCenterOfRotation/1000; % [m]
chord  = Blade_data.ChordLength/1000; % [m]
twist  = deg2rad(Blade_data.BladeTwist); % [rad]
airfoil_type = Blade_data.Airfoil;

% Calculate swept area of turbine blades
area_swept = pi*(R^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining coefficients of power and thrust (deliverable 1)
wind_vel = 10; % [m/s]
rot_velocity = 14*(pi/30); % [rad/s]
pitch = 0; % [deg]

% Calculate angular induction factors and speed ratios using function
[alpha_prime_r, lambda_r, lambda] = calculateInductionFactors(radius, ...
    wind_vel, rot_velocity, alpha);

% Define number of sections in the blade's data
sections = height(Blade_data);

% Calculate drag, lift, z, and theta  forces per unit length
[FL_prime, FD_prime, Fz_prime, Ftheta_prime, Urel, AoA, CD_r, CL_r] = ...
calculateDragLift(wind_vel, rot_velocity, radius, alpha, alpha_prime_r, ...
pitch, twist, chord, dens_air, visc_air, sections, airfoil_type, airfoils);

% Determine force generated from blades
% integration of z force/length along blade
force_1_blade = cumtrapz(radius, Fz_prime); 
% integration of theta force/length along blade multiplied by radius
torque_1_blade = cumtrapz(radius, (radius.*Ftheta_prime));

F = n_blades*max(force_1_blade); % multiplying force per blade by 3 blades
T = n_blades*max(torque_1_blade);% multiplying torque per blade by 3 blades

% Calculating coefficient of thrust and power 
Coeff_thrust = 2*F/(dens_air*(wind_vel^2)*area_swept); 
Coeff_power = 2*T*rot_velocity/(dens_air*(wind_vel^3)*area_swept);

% Report deliverable 1 results 
fprintf('\nDeliverable 1 Results:\n');
fprintf('Coefficent of Power: %.3f\n', Coeff_power);
fprintf('Coefficent of Thrust: %.3f\n', Coeff_thrust);

%% Plots for deliverable 1
% Plotting drag and lift coefficients across blade radius
% figure()
% plot(radius, CD_r, 'k-.', 'LineWidth', 1.5)
% hold on
% plot(radius, CL_r, 'k-','LineWidth', 1.5)
% xlim([0 48])
% xticks(0:4:48)
% xlabel('Blade Radius (m)')
% ylim([0 1.6])
% ylabel('Coeffcient Value')
% legend('Drag Coeffcient', 'Lift Coeffcient', 'Location','east')

% % plotting torque and force across blade radius 
% figure()
% subplot(2,1,1)
% plot(radius, force_1_blade/1000, 'k-','LineWidth', 1.5)
% xlim([0 48])
% xticks(0:4:48)
% xlabel('Blade Radius (m)')
% ylim([0 100])
% ylabel('Force along Blade (kN)')

% subplot(2,1,2)
% plot(radius, torque_1_blade/1000, 'k-','LineWidth', 1.5)
% xlim([0 48])
% xticks(0:4:48)
% xlabel('Blade Radius (m)')
% ylim([0 400])
% ylabel('Torque along Blade (kNm)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimizing blade pitch angle for max CP, with assigned tip speed ratio 
%% and velocity (deliverable 2)

% Deliverable 2 specific parameter
wind_vel_d2 = 7.6; % Assigned wind velocity [m/s]
lambda_d2 = 6.61; % Assigned tip speed ratio 

% Determining the rotational velocity based on tip speed ratio and wind 
% velocity [rad/s]
rot_velocity_d2 = (lambda_d2*wind_vel_d2/R); 

% Varying pitch angles from -30 to 30 degrees and creating empty arrays of
% CP and CT to be analyzed
pitch_values = linspace(-30, 30, 250);
CP_values = zeros(size(pitch_values));
CT_values = zeros(size(pitch_values));

% Calculate induction factors for deliverable 2 conditions
[alpha_prime_r_d2, ~, ~] = calculateInductionFactors(radius, ...
    wind_vel_d2, rot_velocity_d2, alpha);

% Loop to cycle through pitch angles and calculate CP and CT values
for i=1:length(pitch_values)
    pitch_d3 = pitch_values(i); % pitch value that is being looped through

    % Calculate drag, lift, z, and theta  forces per unit length
    [~,~, Fz_prime_d2, Ftheta_prime_d2,~,~,~,~] = calculateDragLift(...
        wind_vel_d2, rot_velocity_d2, radius, alpha, alpha_prime_r_d2, ...
        pitch_d3, twist, chord, dens_air, visc_air, sections, ...
        airfoil_type, airfoils);

    % Determine force generated from blades
    force_1_blade_d2 = cumtrapz(radius, Fz_prime_d2);
    torque_1_blade_d2 = cumtrapz(radius, (radius.*Ftheta_prime_d2));
    % multiplying force and torque per blade by 3 blades
    F_d2 = n_blades*max(force_1_blade_d2); 
    T_d2 = n_blades*max(torque_1_blade_d2);

    % Filling CT and CP arrays with resulting values from pitch angle sweep
    CT_values(i) = 2*F_d2/(dens_air*(wind_vel_d2^2)*area_swept);
    CP_values(i) = 2*T_d2*rot_velocity_d2/(dens_air* ...
        (wind_vel_d2^3)*area_swept);
end

% Determining pitch angle resulting in maximum coeffecient of power
[max_CP, optimized_index] = max(CP_values);
optimal_pitch = pitch_values(optimized_index);
CT_d2 = CT_values(optimized_index);

% Report deliverable 2 results
fprintf(['\nDeliverable 2 Optimization Results (Wind Velocity = ' ...
    '7.6 m/s, Tip Speed Ratio = 6.6):\n']);
fprintf('Optimal Pitch Angle: %.2f degrees\n', optimal_pitch);
fprintf('Maximum Power Coefficient: %.3f\n', max_CP);
fprintf('Resulting Thrust Coefficient: %.3f\n', CT_d2);

% Plot results of optimization
figure()
plot(pitch_values, CP_values, 'k-', 'LineWidth', 1.5)
hold on
plot(optimal_pitch, max_CP, 'ko-','LineWidth', 1, 'MarkerFaceColor','k')
text(optimal_pitch, max_CP, '  Optimized Pitch Angle', ...
    'VerticalAlignment', 'middle')
xlim([-30 30])
ylim([0 0.5])
xticks(-30:5:30)
xlabel('Pitch Angle \beta (degrees)')
ylabel('Coefficient of Power (C_P)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Varying blade pitch angle and tip speed ratio for max CP for an 
%% assigned wind velocity (deliverable 3)

% Varying pitch angle and tip speed ratio variables
wind_vel_d3 = 6.9; % Assigned wind velocity [m/s]
pitch_values_d3 = linspace(-10, 10, 100); % Varying pitch from -10 to 10
lambda_values_d3 = linspace(5, 15, 100); % Varying lambda from 5 to 15

% Generating grid for contour plot
[pitch_grid, lambda_grid] = meshgrid(pitch_values_d3, lambda_values_d3);

% Creating empty arrays of CP and CT to be analyzed
CP_values_d3 = zeros(size(pitch_grid));
CT_values_d3 = zeros(size(pitch_grid));

% Calculating CP for each pitch angle and tip speed for contour plot
for i = 1:height(pitch_grid) % Looping through various tip speed ratios
    lambda_d3 = lambda_grid(i, 1);
    %calculating rotational velocity from tip speed ratio
    omega_d3 = (lambda_d3*wind_vel_d3)/R; 

    % Calculating induction factors for various tip speed ratios
    [alpha_prime_r_d3, ~, ~] = calculateInductionFactors(radius, ...
        wind_vel_d3, omega_d3, alpha);

    for j = 1:width(pitch_grid) % Looping through various pitch angles
        pitch_d3 = pitch_grid(i, j);

        % Calculate drag, lift, z, and theta  forces per unit length
        [~, ~, Fz_prime_d3, Ftheta_prime_d3, ~, ~, ~, ~] = ...
        calculateDragLift(wind_vel_d3, omega_d3, radius, alpha, ...
        alpha_prime_r_d3, pitch_d3, twist, chord, dens_air, visc_air, ...
        sections, airfoil_type, airfoils);

        % Determine force generated from blades
        % integration of z force/length along blade
        force_1_blade_d3 = cumtrapz(radius, Fz_prime_d3); 
        % integration of the force/length along blade multiplied by radius
        torque_1_blade_d3 = cumtrapz(radius, (radius.*Ftheta_prime_d3)); 
        F_d3 = n_blades*max(force_1_blade_d3); % Total turbine force
        T_d3 = n_blades*max(torque_1_blade_d3); % Total turbine torque

        % Filling CT and CP arrays with resulting values 
        CT_values_d3(i,j) = 2*F_d3/(dens_air*(wind_vel_d3^2)*area_swept);
        CP_values_d3(i,j) = 2*T_d3*omega_d3/...
        (dens_air*(wind_vel_d3^3)*area_swept);
    end
end

% Performing optimization to determine optimal pitch angle and tip speed ratio
global x_steps_d3;

% Converting rotor speed limits to tip speed ratio limits
omega_limits = rotor_speed_range*(pi/30); % [rad/s]
lambda_limits = (omega_limits*R)/wind_vel_d3;

% Upper and lower bounds of optimization
lowerbounds_d3 = [-10, lambda_limits(1)];  %[minimum pitch, minimum lambda]
upperbounds_d3 = [ 10, lambda_limits(2)]; %[maximum pitch, maximum lambda]

% Initial guess 
x0_d3 = [0, 8]; % [pitch, lambda]
x_steps_d3 = x0_d3; % Record step iteration history

options_d3 = optimoptions('fmincon', 'Display', 'off', 'Algorithm', ...
    'sqp','OutputFcn', @outfun_d3, 'DiffMinChange', 1e-3);

% Minimizing -CP to maximize CP as we are using 'fmincon' (a minimization
% optimization function). Nagtive CP flips the CP values, so the min is now
% the max
[x_opt_d3, neg_CP_d3, exitflag_d3, output_d3] = fmincon(...
    @(x) obj_neg_CP(x, wind_vel_d3, R, radius, alpha, twist, chord, ...
    dens_air, visc_air, sections, airfoil_type, airfoils, n_blades, ...
    area_swept),x0_d3, [], [], [], [], lowerbounds_d3, upperbounds_d3, ...
    [], options_d3);

% Record optimal pitch angle and tip speed ratio
optimal_pitch_d3  = x_opt_d3(1);
optimal_lambda_d3 = x_opt_d3(2);

% Maximium coefficient of power
max_CP_d3 = -neg_CP_d3;

% Calculating resulting Coefficient of Thrust
% resulting rotational velocity
opt_omega_d3 = (optimal_lambda_d3*wind_vel_d3)/R; 

% Calculating resulting forces on the turbine using optimum values
[opt_ap_d3, ~, ~] = calculateInductionFactors(radius, wind_vel_d3, ...
    opt_omega_d3, alpha);
[~, ~, opt_Fz_d3, ~, ~, ~, ~, ~] = calculateDragLift( ...
    wind_vel_d3, opt_omega_d3, radius, alpha, opt_ap_d3, ...
    optimal_pitch_d3, twist, chord, dens_air, visc_air, sections, ...
    airfoil_type, airfoils);

% resulting force on the turbine
opt_F_d3 = n_blades*max(cumtrapz(radius, opt_Fz_d3)); 
%resulting CT
opt_CT_d3 = (2*opt_F_d3)/(dens_air*(wind_vel_d3^2)*area_swept); 


% Report deliverable 3 results
fprintf(['\nDeliverable 3 Optimization Results (Wind Velocity =' ...
    ' 6.9 m/s):\n']);
fprintf('Iterations: %d\n', output_d3.iterations);
fprintf('Function Evaluations: %d\n', output_d3.funcCount )
fprintf('Optimal Pitch Angle: %.2f degrees\n', optimal_pitch_d3);
fprintf('Optimal Tip Speed Ratio: %.2f\n', optimal_lambda_d3);
fprintf('Maximum Power Coefficient: %.3f\n', max_CP_d3);
fprintf('Resulting Thrust Coefficient: %.3f\n', opt_CT_d3);

% Contour plot over to visualize optimization
figure('Color', 'w');

 % contour with 15 contours
[C_d3, h_d3] = contour(pitch_grid, lambda_grid, CP_values_d3, 15, ...
    'LineWidth', 1.2);
ylabel(colorbar, 'Coefficient of Power (C_P)');
hold on;
clabel(C_d3, h_d3, 'Color', 'k','Fontsize',8); %labeling contours

% Plotting rotor speed constraints
yline(lambda_limits(1), 'm--', 'LineWidth', 1.5); % minimum rotor speed
yline(lambda_limits(2), 'm-', 'LineWidth', 1.5); % maximum rotor speed

% Plot optimization path, initial guess, and optimal CP
% optimziation path
plot(x_steps_d3(:,1), x_steps_d3(:,2), 'r-', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 5, 'LineWidth', 2);
% initial gues
plot(x0_d3(1), x0_d3(2), 'bo', 'MarkerFaceColor', 'b');
% optimal CP
plot(optimal_pitch_d3, optimal_lambda_d3, 'o', 'MarkerFaceColor', 'y', ...
    'MarkerSize', 10);

xlabel('Blade Pitch Angle \beta (degrees)');
ylabel('Tip Speed Ratio \lambda');
xlim([-10 10]);
ylim([5 15]);
grid on;
legend('C_P Contours', 'Min Rotor Speed (9.6 RPM)', ...
    'Max Rotor Speed (15.5 RPM)','Optimization Path', 'Initial Guess', ...
    'Optimal Design Point', 'Location', 'NorthEast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deliverable 4: Determine Required Blade Pitch Angle

% Define the wind speed range from rated to cut-out speed (m/s)
wind_speeds_d4 = 11:0.5:25;

% Define the range of blade pitch angles to test (degrees)
pitch_angles_d4 = 0:0.1:360;

% Set initial values for the required pitch angle and maximum power
required_pitch_d4 = NaN;
max_power_d4 = NaN;

% Test each pitch angle, starting from the smallest
for pitch = pitch_angles_d4

    % Reset maximum power for the current pitch angle
    maximum_power = 0;

    % Check turbine power at each wind speed
    for V = wind_speeds_d4

        % Convert the constant rotor speed of 15.5 RPM to rad/s
        omega = 15.5*pi/30;

        % Calculate turbine power at the current wind speed and pitch
        P = calculateTurbinePower(V, omega, pitch, ...
            radius, chord, twist, alpha, dens_air, ...
            visc_air, sections, airfoil_type, airfoils, n_blades);

        % Stop testing this pitch angle if the power is invalid
        if ~isfinite(P)
            maximum_power = NaN;
            break
        end

        % Update the maximum power for this pitch angle
        maximum_power = max(maximum_power, P);
    end

    % Check whether the maximum power is within the rated limit
    if isfinite(maximum_power) && maximum_power <= power_rated

        % Save the first pitch angle that meets the power limit
        required_pitch_d4 = pitch;
        max_power_d4 = maximum_power;
        break
    end
end

% Display the results
if isnan(required_pitch_d4)
    fprintf('No valid pitch angle found.\n');
else
    fprintf('\nDeliverable 4 Results\n');
    fprintf('Required pitch angle: %.1f degrees\n', ...
        required_pitch_d4);
    fprintf('Maximum power in range: %.4f MW\n', ...
        max_power_d4 / 1e6);
end

%% Deliverable 4 (Group 29 Specific Wind Speed: 16.4 m/s) -> Feeds Deliverable 5
U_hub_d4 = 16.4;                    % [m/s] Group 29 assigned wind speed (Q4 & Q5)
omega_d4 = 15.5 * (pi / 30);        % [rad/s] maximum rotor speed (15.5 RPM)

required_pitch_g29 = NaN;
F_thrust_d4        = NaN;
P_g29              = NaN;

for pitch_test = 0:0.05:45
    [ap_d4, ~, ~] = calculateInductionFactors(radius, U_hub_d4, omega_d4, alpha);
    [~, ~, Fz_d4, Ftheta_d4, ~, ~, ~, ~] = calculateDragLift( ...
        U_hub_d4, omega_d4, radius, alpha, ap_d4, pitch_test, twist, chord, ...
        dens_air, visc_air, sections, airfoil_type, airfoils);

    if any(~isfinite(Ftheta_d4)) || any(~isfinite(Fz_d4))
        continue
    end

    T_total = n_blades * max(cumtrapz(radius, radius .* Ftheta_d4));
    P_test  = T_total * omega_d4;

    if P_test <= power_rated
        required_pitch_g29 = pitch_test;
        P_g29              = P_test;
        F_thrust_d4        = n_blades * max(cumtrapz(radius, Fz_d4)); % [N] Thrust for Deliv 5!
        break
    end
end

fprintf('\nDeliverable 4 Results (Group 29: U = %.1f m/s at 15.5 RPM):\n', U_hub_d4);
fprintf('Required Pitch Angle : %.2f degrees\n', required_pitch_g29);
fprintf('Resulting Power      : %.4f MW\n', P_g29 / 1e6);
fprintf('Total Rotor Thrust   : %.2f kN (Used in Deliverable 5)\n', F_thrust_d4 / 1e3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local Functions
% Calculates local speed ratio, lambda_r, tip speed ratio, lambda, and
% local angular induction factor, alpha_prime
function [alpha_prime_r, lambda_r, lambda] = ...
    calculateInductionFactors(radius, wind_vel, rot_velocity, alpha)
   
    lambda_r = (rot_velocity.*radius)/wind_vel;
    lambda = max(lambda_r);
    
    alpha_prime_r = 0.5*(-1+sqrt(1+(4./(lambda_r.^2))*alpha*(1-alpha)));
end

% Calculates the lift, drag, z, and theta forces per length of the blade
function [FL_prime, FD_prime, Fz_prime, Ftheta_prime, Urel, AoA,... 
    CD_values, CL_values] = calculateDragLift(U, omega, r, alpha, ...
    alpha_prime, pitch, twist, chord, density, viscosity, ...
    sections, names, airfoils)
    
    % Calculating velocities along z and theta axes
    Uz = U.*(1 - alpha); % Velocity along z axis
    Utheta = r.*omega.*(1 + alpha_prime);  % Velocity along theta axis 
    
    % Calculating the magnitude of the relative wind velocity
    Urel = sqrt(Uz.^2 + Utheta.^2);

    % Calculating phi,from inverse tangent of Uz and Utheta
    phi = atan2(Uz, Utheta);

    % Angle of Attack 
    AoA = phi - deg2rad(pitch) - twist; 
    AoA_deg = rad2deg(AoA); % Convert AoA to degrees
    
    % Calculating drag and lift coefficients for each section of the blade
    CL_values=zeros(sections, 1);%Array to fill with drag coeffcient values
    CD_values=zeros(sections, 1);%Array to fill with lift coeffcient values

    % Loop of each section along the blade to select airfoil data and
    % determine coefficients of lift and drag

    for i = 1:sections
        section_airfoil = char(names(i));
        c_AoA = AoA_deg(i); % AoA that is being looped through
    
        % Determine CD and CL for cylindrical blade sections using
        % Reynold's Number 
        if contains(section_airfoil, 'circle')
            CL_values(i) = 0; % Lift is zero at Cylindrical section
            
            % Compute Reynolds number based on relative velocity and
            % characteristic legnth (chord distance)
            Re = (density*Urel(i)*chord(i))/viscosity;
            
            % Calculate Drag Coefficient via drag coefficient function
            CD_values(i) = cylinderCD(Re);

        % Matching airfoil tables and storing CL and CD
        elseif contains(section_airfoil, 'DU 91-W2-250')
            airfoil_table = airfoils.DU91;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, ...
                c_AoA, 'linear','extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, ...
                c_AoA, 'linear','extrap');
            
        elseif contains(section_airfoil, 'DU 93-W-210')
            airfoil_table = airfoils.DU93;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, ...
                c_AoA, 'linear','extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, ...
                c_AoA, 'linear','extrap');
            
        elseif contains(section_airfoil, 'DU 96-W-180')
            airfoil_table = airfoils.DU96;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, ...
                c_AoA, 'linear','extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, ...
                c_AoA, 'linear','extrap');
            
        elseif contains(section_airfoil, 'DU 97-W-300')
            airfoil_table = airfoils.DU97;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, ...
                c_AoA, 'linear','extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, ...
                c_AoA, 'linear','extrap');
        end
    end
    
    % Determing drag and lift forces/length along blade
    % Lift force/length [N/m]
    FL_prime = 0.5*density.*(Urel.^2).*chord.*CL_values;
    % Drag force/length [N/m]
    FD_prime = 0.5*density.*(Urel.^2).*chord.*CD_values;
    
    % Converting Lift and Drag into z and theta forces 
    Fz_prime = FL_prime.*cos(phi) + FD_prime.*sin(phi); 
    Ftheta_prime = FL_prime.*sin(phi) - FD_prime.*cos(phi);
end


function [C_D] = cylinderCD(Re)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: cylinderCD
%
% PURPOSE:
% To calculate and return the coeffecient of drag for a cyliner in cross
% flow.
%
% Fit functions can be found in references:
%
% - White, F.M. (2006) Viscous Fluid Flow. 3rd Edition, McGraw-Hill, Boston.
% - Nian-Sheng Cheng, *Calculation of Drag Coefficient for Arrays of ...
% Emergent Circular Cylinders with Pseudofluid Model.*...
% J. Hydraul. Eng. 2013.139:602-611.
%
% INPUT
% Re - Reynolds number [-]
%
% OUTPUT
% C_D - drag coeffecient [-]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: MJA
% DATE: 2024.10.21
%
% DESCRIPTION OF LOCAL VARIABLES
%
%
% FUNCTIONS CALLED
% log10 (MATLAB)
% tanh (MATLAB)
%
% START OF EXECUTABLE CODE
%
% piecewise fit functions based on various reynolds number regimes for a
% cylinder in cross-flow
    if Re < 2*10^5
        C_D = 11 * Re.^(-0.75) + 0.9 * (1.0 - exp(-1000./Re))...
        + 1.2 * (1.0 -exp(-(Re./4500).^0.7));
    elseif Re <= 5*10^5
        C_D = 10.^(0.32*tanh(44.4504 - 8 * log10 (Re))-0.238793158);
    else
        C_D = 0.1 * log10(Re) - 0.2533429;
    end
end

% Function to record iteration steps
function stop = outfun_d3(x, optimValues, state)
    stop = false; % Keeps solver running
    global x_steps_d3;
    if strcmp(state, 'iter')
        % Stores pitch angle and tip speed ratio at each step
        x_steps_d3 = [x_steps_d3; x]; 
    end
end

% Function to calculate the negative CP values for deliverable 3, minimal
% negative CP is maximum CP
function [neg_CP, opt_CT] = obj_neg_CP(x, wind_vel, R, radius, ...
    alpha, twist, chord, dens_air, visc_air, sections, airfoil_type, ...
    airfoils, n_blades, area_swept)
    % pitch angle, tip speed ratio, and rotational velocity values
    pitch = x(1);
    lambda = x(2);
    omega = (lambda*wind_vel)/R;
    
    % Calculating angle induction factors and forces on the turbine blade
    % using functions
    [alpha_prime_r, ~, ~] = calculateInductionFactors(radius, wind_vel, ...
        omega, alpha);
    [~, ~, Fz_prime, Ftheta_prime, ~, ~, ~, ~] = calculateDragLift( ...
    wind_vel, omega, radius, alpha, alpha_prime_r, pitch, twist, chord, ...
    dens_air, visc_air, sections, airfoil_type, airfoils);
        
    % Determining coefficient of power
    force_1_blade = cumtrapz(radius, Fz_prime); % force in 1 blade
    torque_1_blade =cumtrapz(radius,(radius.*Ftheta_prime));%torque 1 blade
    F = n_blades*max(force_1_blade); %force over the whole turbine
    T = n_blades*max(torque_1_blade); %torque over the whole turbine
    CP = (2*T*omega)/(dens_air*(wind_vel^3)*area_swept); %CP
    CT = 2*F/(dens_air*(wind_vel^2)*area_swept); % CT
    
    neg_CP = -CP; % Returning negative CP to maximize CP
    opt_CT = CT; %Returning resulting CT
    
end

function P = calculateTurbinePower(V, omega, pitch, ...
    radius, chord, twist, alpha, density, viscosity, ...
    sections, airfoil_type, airfoils, n_blades)

    % Calculate the induction factors along the blade
    [alpha_prime_r, ~, ~] = calculateInductionFactors( ...
        radius, V, omega, alpha);

    % Calculate the tangential force per unit blade length
    [~, ~, ~, Ftheta_prime, ~, ~, ~, ~] = ...
        calculateDragLift(V, omega, radius, alpha, ...
        alpha_prime_r, pitch, twist, chord, density, ...
        viscosity, sections, airfoil_type, airfoils);

    % Return NaN if any tangential force value is invalid
    if any(~isfinite(Ftheta_prime))
        P = NaN;
        return
    end

    % Integrate the tangential force along the blade to find torque
    torque_1_blade = cumtrapz(radius, radius .* Ftheta_prime);

    % Calculate total rotor torque from all blades
    torque = n_blades * max(torque_1_blade);

    % Calculate power from torque and angular velocity (W)
    P = torque * omega;
end