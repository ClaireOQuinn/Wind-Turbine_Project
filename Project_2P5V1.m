%% Wind Turbine Project
clc
clear
close all

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

%% Determining coefficient of power and coefficient of thrust (deliverable 1)
wind_vel = 10; % [m/s]
rot_velocity = 14*(pi/30); % [rad/s]
pitch = 0; % [deg]

% Calculate angular induction factors and speed ratios using function
[alpha_prime_r, lambda_r, lambda] = calculateInductionFactors(radius, wind_vel, rot_velocity, alpha);

% Define number of sections in the blade's data
sections = height(Blade_data);

% Calculate drag, lift, z, and theta forces per unit length
[FL_prime, FD_prime, Fz_prime, Ftheta_prime, Urel, AoA, CD_r, CL_r] = calculateDragLift( ...
    wind_vel, rot_velocity, radius, alpha, alpha_prime_r, pitch, twist, chord, ...
    dens_air, visc_air, sections, airfoil_type, airfoils);

% Determine force generated from blades
force_1_blade = cumtrapz(radius, Fz_prime); % integration of z force/length along blade
torque_1_blade = cumtrapz(radius, (radius.*Ftheta_prime)); % integration of theta force/length along blade multiplied by radius

F = n_blades*max(force_1_blade); % multiplying force per blade by 3 blades
T = n_blades*max(torque_1_blade);  % multiplying torque per blade by 3 blades

% Calculating coefficient of thrust and power 
Coeff_thrust = 2*F/(dens_air*(wind_vel^2)*area_swept); % Coefficient of Thrust
Coeff_power = 2*T*rot_velocity/(dens_air*(wind_vel^3)*area_swept); % Coefficient of Power

% Report deliverable 1 results 
fprintf('\nDeliverable 1 Optimization Results:\n');
fprintf('Coefficent of Power: %.3f\n', Coeff_power);
fprintf('Coefficent of Thrust: %.3f\n', Coeff_thrust);

%% Plots for deliverable 1
% Plotting drag and lift coefficients across blade radius
figure()
plot(radius, CD_r, 'k-.', 'LineWidth', 1.5)
hold on
plot(radius, CL_r, 'k-','LineWidth', 1.5)
xlim([0 48])
xticks(0:4:48)
xlabel('Blade Radius (m)')
ylim([0 1.6])
ylabel('Coeffcient Value')
legend('Drag Coeffcient', 'Lift Coeffcient', 'Location','east')

% plotting torque and force across blade radius 
figure()
subplot(2,1,1)
plot(radius, force_1_blade/1000, 'k-','LineWidth', 1.5)
xlim([0 48])
xticks(0:4:48)
xlabel('Blade Radius (m)')
ylim([0 100])
ylabel('Force along Blade (kN)')

subplot(2,1,2)
plot(radius, torque_1_blade/1000, 'k-','LineWidth', 1.5)
xlim([0 48])
xticks(0:4:48)
xlabel('Blade Radius (m)')
ylim([0 400])
ylabel('Torque along Blade (kNm)')

%% Optimizing blade pitch angle for max CP, with assigned tip speed ratio and velocity (deliverable 2)

% Deliverable 2 specific parameter
wind_vel_d2 = 7.6; % Assigned wind velocity [m/s]
lambda_d2 = 6.61; % Assigned tip speed ratio 
rot_velocity_d2 = (lambda_d2*wind_vel_d2/R); % Determining the rotational velocity based on tip speed ratio and wind velocity [rad/s]

% Varying pitch angles from -30 to 30 degrees and creating empty arrays of
% CP and CT to be analyzed
pitch_values = linspace(-30, 30, 250);
CP_values = zeros(size(pitch_values));
CT_values = zeros(size(pitch_values));

% Calculate induction factors for deliverable 2 conditions
[alpha_prime_r_d2, ~, ~] = calculateInductionFactors(radius, wind_vel_d2, rot_velocity_d2, alpha);

% Loop to cycle through pitch angles and calculate CP and CT values
for j=1:length(pitch_values)
    c_pitch = pitch_values(j); % pitch value that is being looped through

    % Calculate drag, lift, z, and theta  forces per unit length
    [~, ~, Fz_prime_d2, Ftheta_prime_d2, ~, ~, ~, ~] = calculateDragLift( ...
        wind_vel_d2, rot_velocity_d2, radius, alpha, alpha_prime_r_d2, ...
        c_pitch, twist, chord, dens_air, visc_air, sections, airfoil_type, airfoils);

    % Determine force generated from blades
    force_1_blade_d2 = cumtrapz(radius, Fz_prime_d2); % integration of z force/length along blade
    torque_1_blade_d2 = cumtrapz(radius, (radius.*Ftheta_prime_d2)); % integration of theta force/length along blade multiplied by radius
    F_d2 = n_blades*max(force_1_blade_d2); % multiplying force per blade by 3 blades
    T_d2 = n_blades*max(torque_1_blade_d2); % multiplying force per blade by 3 blades

    % Filling CT and CP arrays with resulting values from pitch angle sweep
    CT_values(j) = 2*F_d2/(dens_air*(wind_vel_d2^2)*area_swept);
    CP_values(j) = 2*T_d2*rot_velocity_d2/(dens_air*(wind_vel_d2^3)*area_swept);
end

% Determining pitch angle resulting in maximum coeffecient of power
[max_CP, optimized_index] = max(CP_values);
optimal_pitch = pitch_values(optimized_index);
CT_d2 = CT_values(optimized_index);

% Report deliverable 2 results
fprintf('\nDeliverable 2 Optimization Results:\n');
fprintf('Optimal Pitch Angle: %.2f deg\n', optimal_pitch);
fprintf('Maximum Power Coefficient: %.3f\n', max_CP);
fprintf('Resulting Thrust Coefficient: %.3f\n', CT_d2);

% Plot results of optimization
figure()
plot(pitch_values, CP_values, 'k-', 'LineWidth', 1.5)
hold on
plot(optimal_pitch, max_CP, 'ko-','LineWidth', 1, 'MarkerFaceColor','k')
text(optimal_pitch, max_CP, '  Optimized Pitch Angle', 'VerticalAlignment', 'middle')
xlim([-30 30])
ylim([0 0.5])
xticks(-30:5:30)
xlabel('Pitch Angle (degrees)')
ylabel('Coefficient of Power')

%% Varying blade pitch angle and tip speed ratio for max CP for an assigned wind velocity (deliverable 3)

% Varying pitch angles from -30 to 30 and varying lambda from 1 to 10,
wind_vel_d3 = 6.9; % Assigned wind velocity [m/s]
pitch_values_d3 = linspace(-30, 30, 250); 
lambda_values_d3 = linspace(1, 10, 250);

% Creating empty arrays of CP and CT to be analyzed
CP_values_d3 = zeros(size(pitch_values_d3));
CT_values_d3 = zeros(size(pitch_values_d3));

%% =========================================================================
%% DELIVERABLE 5 - PART 1: TOWER DEFORMATION & STATIC BENDING YIELD (ODE45)
%% =========================================================================
% Euler-Bernoulli Beam Assumptions:
% 1. Pure bending in the x-y plane (front-back / up-down; no side-to-side)
% 2. Transverse shear strain = 0 (plane sections remain plane)
% 3. Narrow beam approximation (transverse/hoop stress negligible)
% 4. Free outer surface (zero radial/surface stress on outer tower wall)

% --- 1. Structural Constants & Component Masses ---
g = 9.81; % gravity [m/s^2]
E_steel   = 200e9;    % [Pa] Young's modulus for ASTM A572 Grade 50 steel
m_rotor   = 55.7e3;   % [kg] rotor assembly mass 
m_blade   = 10.7e3;   % [kg] single blade mass
m_nacelle = 94.3e3;   % [kg] nacelle mass
m_top     = m_rotor + m_nacelle + 3*m_blade; % [kg] total Rotor-Nacelle Assembly (RNA) mass

% --- 2. Deliverable 4 Operating Conditions (Group 29) ---
U_hub_d4  = 16.4;     % [m/s] Group 29 assigned wind speed for Q4/Q5
alpha_abl = 1/7;      % [-] atmospheric boundary layer exponent (Assumption 4)
F_thrust  = 350e3;    % [N] <-- REPLACE with team's Deliverable 4 total thrust F_d4

% --- 3. Run Static Deflection (ODE45) & Bending Stress Analysis ---
[z_tower, w_static, w_hub_static, sigma_axial, sigma_bending, sigma_static_max, ...
 SF_yield_static, M_bend, V_shear, P_axial, A_tower, I_tower, OD_tower, m_tower_total] = ...
    analyzeTowerStatic(tower_data, hub_height, U_hub_d4, F_thrust, m_top, ...
    dens_air, visc_air, dens_steel, E_steel, strength_yield_steel, alpha_abl, g);

% Display Static Results
[peak_static_stress, idx_peak] = max(sigma_static_max);
fprintf('\n--- Deliverable 5 (Part 1): Static Deflection & Stress (ODE45) ---\n');
fprintf('Total Tower Mass        : %.1f metric tons\n', m_tower_total / 1000);
fprintf('Base Bending Moment     : %.2f MN*m\n', M_bend(1) / 1e6);
fprintf('Max Static Normal Stress: %.2f MPa at z = %.2f m\n', peak_static_stress / 1e6, z_tower(idx_peak));
fprintf('Static Yield Safety Fac.: %.2f (Yield = 345 MPa)\n', SF_yield_static);
fprintf('Tower Top Deflection    : %.3f m (%.1f mm) at z = 77.7 m\n', w_static(end), w_static(end) * 1000);
fprintf('Hub Center Deflection   : %.3f m (%.1f mm) at z = 80.4 m\n', w_hub_static, w_hub_static * 1000);

% --- 4. Required Report Plots: Static Deflection & Bending Stress ---
% Plot 1: Static Deflection Plot (Q5) - Extended with rigid nacelle link to 80.4 m
z_full = [z_tower; hub_height];
w_full = [w_static; w_hub_static];

figure('Name', 'Q5: Static Deflection Plot', 'Color', 'w')
plot(w_full * 1000, z_full, 'b-', 'LineWidth', 1.8)
hold on
yline(z_tower(end), 'k--', 'Tower Top Flange (77.7 m)', 'LabelHorizontalAlignment', 'left')
xlabel('Lateral Deflection (mm)')
ylabel('Height z (m)')
title('Static Tower Deflection vs. Height (ODE45)')
grid on

% Plot 2: Bending Stress - Tower Height Plot (Q5)
figure('Name', 'Q5: Bending Stress - Tower Height Plot', 'Color', 'w')
plot(sigma_bending / 1e6, z_tower, 'r-', 'LineWidth', 1.8)
hold on
xline(strength_yield_steel / 1e6, 'k--', 'Yield Strength (345 MPa)', 'LineWidth', 1.5)
xlabel('Bending Stress (MPa)')
ylabel('Tower Height z (m)')
title('Static Bending Stress vs. Tower Height')
legend('Bending Stress', 'Yield Strength (345 MPa)', 'Location', 'northeast')
grid on

%% =========================================================================
%% DELIVERABLE 5 - PART 2: 1-DOF TOWER VIBRATION & STRESS SUPERPOSITION
%% =========================================================================
omega_d4   = rotor_speed_range(2) * (pi / 30); % [rad/s] max rotor speed (15.5 RPM)
zeta_tower = 0.015;                            % [-] structural damping ratio (1.5% for steel)
nHarmonics = 50;                               % [-] number of Fourier harmonics (< 500)

[t_vib, x_dyn, F_shadow_pulse, k_tower, m_eq, omega_n, f_n, ...
 x_dyn_amp, w_total_max, sigma_total_max, SF_yield_total] = ...
    analyzeTowerVibration(z_tower, OD_tower, I_tower, w_static, sigma_axial, ...
    sigma_bending, R, omega_d4, n_blades, F_thrust, m_top, m_tower_total, ...
    E_steel, strength_yield_steel, zeta_tower, nHarmonics);

fprintf('\n--- Deliverable 5 (Part 2): 1-DOF Vibration & Superposition ---\n');
fprintf('Equivalent Stiffness k  : %.3f MN/m\n', k_tower / 1e6);
fprintf('Natural Frequency wn    : %.3f rad/s (%.3f Hz)\n', omega_n, f_n);
fprintf('3P Blade-Pass Freq.     : %.3f Hz\n', n_blades * rotor_speed_range(2) / 60);
fprintf('Dynamic Deflection Amp. : %.2f mm (Oscillatory)\n', x_dyn_amp * 1000);
fprintf('Peak Dynamic Deflection : %.2f mm (Including shadow DC shift)\n', max(abs(x_dyn)) * 1000);
fprintf('Total Tip Deflection    : %.2f mm (Static + Dynamic)\n', w_total_max * 1000);
fprintf('Total Max Bending Stress: %.2f MPa (Superposition at base)\n', sigma_total_max / 1e6);
fprintf('Total Yield Safety Fac. : %.2f (Yield = 345 MPa)\n', SF_yield_total);

% Plot 3: Dynamic Deflection Plot (Q5 - Required by Rubric)
figure('Name', 'Q5: Dynamic Deflection Plot', 'Color', 'w')
plot(t_vib, x_dyn * 1000, 'b-', 'LineWidth', 1.6)
xlabel('Time t (s)')
ylabel('Dynamic Tip Deflection x_{dyn} (mm)')
title('1-DOF Tower Dynamic Deflection Under Blade Shadow Effect')
grid on

%% ========================================================================
%% Functions (MUST remain at the very bottom of the MATLAB script!)
%% ========================================================================

% Calculates local speed ratio, lambda_r, tip speed ratio, lambda, and
% local angular induction factor, alpha_prime
function [alpha_prime_r, lambda_r, lambda] = calculateInductionFactors(radius, wind_vel, rot_velocity, alpha)
   
    lambda_r = (rot_velocity.*radius)/wind_vel;
    lambda = max(lambda_r);
    
    alpha_prime_r = 0.5*(-1 + sqrt(1 + (4./(lambda_r.^2))*alpha*(1 - alpha)));
end

% Calculates the lift, drag, z, and theta forces per unit length of the blade
function [FL_prime, FD_prime, Fz_prime, Ftheta_prime, Urel, AoA, CD_values, CL_values] = calculateDragLift( ...
    U, omega, r, alpha, alpha_prime, pitch, twist, chord, density, viscosity, sections, names, airfoils)
    
    % Calculating velocities along z and theta axes
    Uz = U.*(1 - alpha); % Velocity along z axis
    Utheta = r.*omega.*(1 + alpha_prime);  % Velocity along theta axis 
    
    % Calculating the magnitude of the relative wind velocity
    Urel = sqrt(Uz.^2 + Utheta.^2); % Pythagoren theorem to calculate the relative velocity magnitude

    % Calculating phi, angle of relative wind
    phi = atan2(Uz, Utheta); % Finding phi from inverse tangent of Uz and Utheta

    % Angle of Attack 
    AoA = phi - deg2rad(pitch) - twist; % Determining AoA (twist is already in radians)
    AoA_deg = rad2deg(AoA); % Convert AoA to degrees
    
    % Calculating drag and lift coefficients for each section of the blade
    CL_values = zeros(sections, 1); % Empty array to fill with drag coefficeint vlaues
    CD_values = zeros(sections, 1); % Empty array to fill with lift coefficeint vlaues

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
            
            % Calculate Drag Coefficient via cylinder drag coefficient function
            CD_values(i) = cylinderCD(Re);

        % Matching airfoil tables and storing CL and CD
        elseif contains(section_airfoil, 'DU 91-W2-250')
            airfoil_table = airfoils.DU91;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear');
            
        elseif contains(section_airfoil, 'DU 93-W-210')
            airfoil_table = airfoils.DU93;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear');
            
        elseif contains(section_airfoil, 'DU 96-W-180')
            airfoil_table = airfoils.DU96;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear');
            
        elseif contains(section_airfoil, 'DU 97-W-300')
            airfoil_table = airfoils.DU97;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear');
        end
    end
    
    % Determing drag force/length along blade and lift force/legnth along
    % blade
    FL_prime = 0.5*density.*(Urel.^2).*chord.*CL_values; % Lift force/length [N/m]
    FD_prime = 0.5*density.*(Urel.^2).*chord.*CD_values; % Drag force/length [N/m]
    
    % Converting Lift and Drag into z and theta forces 
    Fz_prime = FL_prime.*cos(phi) + FD_prime.*sin(phi); % Using trig to determine z force 
    Ftheta_prime = FL_prime.*sin(phi) - FD_prime.*cos(phi); % Using trig to determine theta force
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
    if Re <= 0
        C_D = 0;
    elseif Re < 2*10^5
        C_D = 11 * Re.^(-0.75) + 0.9 * (1.0 - exp(-1000./Re))...
        + 1.2 * (1.0 -exp(-(Re./4500).^0.7));
    elseif Re <= 5*10^5
        C_D = 10.^(0.32*tanh(44.4504 - 8 * log10 (Re))-0.238793158);
    else
        C_D = 0.1 * log10(Re) - 0.2533429;
    end
end


function [z, deflection, w_hub, sigma_axial, sigma_bending, sigma_comp_max, ...
          SF_yield, M_bend, V_shear, P_axial, A, I, OD, m_tower_total] = ...
          analyzeTowerStatic(tower_data, hub_height, U_hub, F_thrust, m_top, ...
          dens_air, visc_air, dens_steel, E_steel, sigma_y, alpha_abl, g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: analyzeTowerStatic
%
% PURPOSE:
% Computes internal shear, bending moment, axial loads, normal stresses,
% and static lateral deflection (via ODE45) for a tapered wind turbine
% tower subjected to rotor thrust and atmospheric boundary layer wind drag.
%
% INPUT
% tower_data - table containing Height (mm), OD (mm), and Wall thk (mm)
% hub_height - rotor hub elevation [m]
% U_hub      - reference wind velocity at hub height [m/s]
% F_thrust   - total aerodynamic rotor thrust at hub [N]
% m_top      - combined rotor and nacelle mass at tower top [kg]
% dens_air   - air density [kg/m^3]
% visc_air   - dynamic viscosity of air [N*s/m^2]
% dens_steel - structural steel density [kg/m^3]
% E_steel    - Young's modulus of steel [Pa]
% sigma_y    - yield strength of steel [Pa]
% alpha_abl  - atmospheric boundary layer power-law exponent [-]
% g          - gravitational acceleration [m/s^2]
%
% OUTPUT
% z              - discretized tower elevation vector [m]
% deflection     - lateral static deflection along tower [m]
% w_hub          - lateral static deflection at hub center [m]
% sigma_axial    - compressive axial stress profile [Pa]
% sigma_bending  - bending normal stress profile at outer fiber [Pa]
% sigma_comp_max - combined peak compressive stress profile [Pa]
% SF_yield       - minimum static yield factor of safety [-]
% M_bend         - internal bending moment profile [N*m]
% V_shear        - internal shear force profile [N]
% P_axial        - internal axial compressive force profile [N]
% A              - cross-sectional steel area profile [m^2]
% I              - area moment of inertia profile [m^4]
% OD             - outer diameter profile [m]
% m_tower_total  - total structural mass of steel tower [kg]
%
% FUNCTIONS CALLED
% interp1, cumtrapz, ode45, trapz, cylinderCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1. Extract and convert raw tower dimensions (mm -> m)
    z_raw = tower_data{:, 1} / 1000; % height of tower
    OD_raw = tower_data{:, 2} / 1000; % OD of tower
    t_raw = tower_data{:, 3} / 1000; % wall thickness

    L_tower  = max(z_raw); % tower height (77.7 m)
    z_offset = hub_height - L_tower; % [m] nacelle offset length (2.7 m)

    % 2. Discretize tower geometry
    N_nodes = 1000;
    z = linspace(0, L_tower, N_nodes)';
    OD = interp1(z_raw, OD_raw, z, 'linear');
    t  = interp1(z_raw, t_raw,  z, 'linear');
    ID = OD - 2 .* t;

    A = (pi / 4)  .* (OD.^2 - ID.^2);
    I = (pi / 64) .* (OD.^4 - ID.^4);

    mu_tower      = dens_steel .* A;
    m_tower_total = trapz(z, mu_tower);

    % 3. Atmospheric Boundary Layer wind profile & distributed tower drag
    U_tower  = U_hub .* (z ./ hub_height).^alpha_abl;
    Re_tower = (dens_air .* U_tower .* OD) ./ visc_air;
    CD_tower = arrayfun(@cylinderCD, Re_tower);
    q_wind   = 0.5 * dens_air .* (U_tower.^2) .* CD_tower .* OD;

    % 4. Integrate internal shear V(z), moment M(z), and axial force P(z)
    V_shear = F_thrust + flipud(cumtrapz(flipud(z), flipud(q_wind))) * (-1);
    M_bend  = (F_thrust * z_offset) + flipud(cumtrapz(flipud(z), flipud(V_shear))) * (-1);

    w_tower = mu_tower .* g;
    P_axial = (m_top * g) + flipud(cumtrapz(flipud(z), flipud(w_tower))) * (-1);

    % 5. Solve Euler-Bernoulli beam deflection ODE via ODE45
    % State vector Y = [deflection w(z); slope dw/dz]
    beamODE = @(z_eval, Y) [ ...
        Y(2); ...
        interp1(z, M_bend, z_eval, 'linear') / (E_steel * interp1(z, I, z_eval, 'linear')) ...
    ];

    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, Y_sol] = ode45(beamODE, z, [0; 0], opts);

    deflection = Y_sol(:, 1);
    slope      = Y_sol(:, 2);
    w_hub      = deflection(end) + slope(end) * z_offset;

    % 6. Compute static normal stresses & yield factor of safety
    sigma_axial    = P_axial ./ A;
    sigma_bending  = (M_bend .* (OD / 2)) ./ I;
    sigma_comp_max = sigma_axial + sigma_bending;
    SF_yield       = sigma_y / max(sigma_comp_max);
end

function [t, x_dyn, F_pulse, k_eq, m_eq, omega_n, f_n, x_osc_amp, ...
          w_total_max, sigma_total_max, SF_yield_total] = ...
          analyzeTowerVibration(z, OD, I, w_static, sigma_axial, sigma_bending, ...
          R, omega_rotor, n_blades, F_thrust, m_top, m_tower, E_steel, sigma_y, zeta, nHarmonics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: analyzeTowerVibration
%
% PURPOSE:
% Models the wind turbine tower as a 1-DOF damped harmonic oscillator
% subjected to periodic blade-passing tower shadow pulses using pulseFFT,
% and computes total combined bending stress via tip deflection superposition.
%
% INPUT
% z, OD, I      - tower height [m], outer diameter [m], moment of inertia [m^4]
% w_static      - static lateral deflection profile [m]
% sigma_axial   - static axial stress profile [Pa]
% sigma_bending - static bending stress profile [Pa]
% R             - rotor blade radius [m]
% omega_rotor   - rotor angular velocity [rad/s]
% n_blades      - number of rotor blades [-]
% F_thrust      - total aerodynamic rotor thrust [N]
% m_top         - mass at tower top (RNA) [kg]
% m_tower       - total mass of structural steel tower [kg]
% E_steel       - Young's modulus of steel [Pa]
% sigma_y       - yield strength of steel [Pa]
% zeta          - structural damping ratio [-]
% nHarmonics    - number of Fourier terms for pulseFFT (< 500)
%
% OUTPUT
% t               - time vector over 4 blade-passing periods [s]
% x_dyn           - dynamic tip deflection time history [m]
% F_pulse         - reconstructed shadow force fluctuation [N]
% k_eq            - equivalent lateral stiffness 3*E*I_avg / L^3 [N/m]
% m_eq            - equivalent modal mass m_top + 0.23*m_tower [kg]
% omega_n, f_n    - natural frequency [rad/s] and [Hz]
% x_osc_amp       - oscillatory dynamic amplitude [m]
% w_total_max     - superimposed static + dynamic max tip deflection [m]
% sigma_total_max - total max normal stress at base [Pa]
% SF_yield_total  - combined yield factor of safety [-]
%
% FUNCTIONS CALLED
% pulseFFT, mean, linspace, cos, sqrt, atan2, max, abs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L_tower = z(end);                   % [m] structural tower height (77.7 m)
    I_avg   = mean(I);                  % [m^4] constant cross-section average I
    D_avg   = mean(OD);                 % [m] constant cross-section average OD

    % 1. Equivalent Stiffness, Rayleigh Mass, and Natural Frequency
    k_eq    = (3 * E_steel * I_avg) / (L_tower^3); % [N/m]
    m_eq    = m_top + 0.23 * m_tower;              % [kg]
    omega_n = sqrt(k_eq / m_eq);                   % [rad/s]
    f_n     = omega_n / (2 * pi);                  % [Hz]

    % 2. Tower Shadow Square-Wave Pulse Parameters
    f_1P        = omega_rotor / (2 * pi);          % [Hz] rotor shaft frequency
    freq_3P     = n_blades * f_1P;                 % [Hz] blade-passing frequency
    T_3P        = 1 / freq_3P;                     % [s] pulse period

    r_eff       = R / 2;                           % [m] mid-span effective radius
    pulseWidth  = (D_avg / r_eff) / omega_rotor;   % [s] time blade shadows tower
    pulseHeight = -F_thrust / n_blades;            % [N] binary loss of 1 blade thrust

    % 3. Fourier Decomposition via pulseFFT
    [fourierMag, fourierPhase] = pulseFFT(freq_3P, pulseWidth, pulseHeight, nHarmonics);

    % 4. Steady-State 1-DOF Harmonic Response
    t       = linspace(0, 4 * T_3P, 2000)';
    x_dc    = fourierMag(1) / k_eq;                % [m] DC mean deflection shift
    x_osc   = zeros(size(t));
    F_pulse = fourierMag(1) * ones(size(t));

    for idx = 2:nHarmonics
        harm_order = idx - 1;                      % harmonic multiple (1X, 2X, ...)
        w_i        = 2 * pi * harm_order * freq_3P;
        r_i        = w_i / omega_n;

        % Reconstruct forcing harmonic
        F_pulse = F_pulse + fourierMag(idx) * cos(w_i * t + fourierPhase(idx));

        % 1-DOF frequency response magnitude & phase lag
        H_mag = (1 / k_eq) / sqrt((1 - r_i^2)^2 + (2 * zeta * r_i)^2);
        psi_i = atan2(2 * zeta * r_i, 1 - r_i^2);

        x_osc = x_osc + (fourierMag(idx) * H_mag) * cos(w_i * t + fourierPhase(idx) - psi_i);
    end

    x_dyn     = x_dc + x_osc;
    x_osc_amp = max(abs(x_osc));

    % 5. Superposition of Static + Dynamic Tip Deflection for Max Stress
    w_stat_tip      = w_static(end);
    x_dyn_peak      = max(abs(x_dyn));
    w_total_max     = w_stat_tip + x_dyn_peak;

    dyn_amp_factor  = 1 + (x_dyn_peak / w_stat_tip);
    sigma_total_max = sigma_axial(1) + sigma_bending(1) * dyn_amp_factor;
    SF_yield_total  = sigma_y / sigma_total_max;
end