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
visc_air = 8e-5; % [N s/m^2] dynamic viscosity of air
dens_steel = 7850; % [kg/m^3] steel density
strength_tens_steel = 450e6; % [Pa] tensile strength of steel
strength_yield_steel = 345e6; % [Pa] yield strength of steel

% Import Turbine Data
Blade_data = readtable('BladeProfile.csv'); % Blade data
tower_data = readtable('towerSpecs.csv'); % Tower data

% Airfoil Data stored in structure 
airfoils = struct();
airfoils.DU91 = readtable('DU91-W2-250.csv');
airfoils.DU93 = readtable('DU93-W-210.csv');
airfoils.DU96 = readtable('DU96-W-180.csv');
airfoils.DU97 = readtable('DU97-W-300.csv');

% Convert mm to m and store blade data 
radius = Blade_data.DistanceFromCenterOfRotation / 1000; % [m]
chord  = Blade_data.ChordLength / 1000;  % [m]
twist  = deg2rad(Blade_data.BladeTwist); % [rad]
airfoil_type = Blade_data.Airfoil;

% Calculate swept area of turbine blades
area_swept = pi*(R^2);

%% Determining coefficient of power and coefficient of thrust (deliverable 1)
wind_vel = 10; % [m/s]
rot_velocity = 14 * (pi/30); % [rad/s]
pitch = 0; % [deg]

% Calculate angular induction factors and speed ratios using function
[alpha_prime_r, lambda_r, lambda] = calculateInductionFactors(radius, wind_vel, rot_velocity, alpha);

% Define number of sections in the blade's data
sections = height(Blade_data);

% Calculate drag, lift, z, and theta  forces per unit length
[FL_prime, FD_prime, Fz_prime, Ftheta_prime, Urel, AoA, CD_r, CL_r] = calculateDragLift( ...
    wind_vel, rot_velocity, radius, alpha, alpha_prime_r, pitch, twist, chord, ...
    dens_air, visc_air, sections, airfoil_type, airfoils);

% Determine force generated from blades
Force_total_blades = n_blades * trapz(radius, Fz_prime); % integration of z force/length along blade, factored by 3 blades
Torque_total_blades = n_blades * trapz(radius, (radius.*Ftheta_prime)); % integration of theta force/length along blade multiplied by radius, factored by 3 blades

Coeff_thrust = 2*Force_total_blades/(dens_air*(wind_vel^2)*area_swept);
Coeff_power = 2*Torque_total_blades*rot_velocity/(dens_air*(wind_vel^2)*area_swept);
fprintf('Coefficent of Power: %.3f\n', Coeff_power);
fprintf('Coefficent of Thrust: %.3f\n', Coeff_thrust);

% Plots for deliverable 1
figure()
plot(radius, CD_r, 'k-.', 'LineWidth', 1.5)
hold on
plot(radius, CL_r, 'k--','LineWidth', 1.5)
xlim([0 48])
xticks([0 4 8 12 16 20 24 28 32 36 40 44 48])
xlabel('Blade Radius (m)')
ylim([0 1.6])
ylabel('Coeffcient Value')
legend('Drag Coeffcient', 'Lift Coeffcient', 'Location','east')
 

%% Functions
function [alpha_prime_r, lambda_r, lambda] = calculateInductionFactors(radius, wind_vel, rot_velocity, alpha)
    % Calculates local speed ratio (lambda_r), tip speed ratio (lambda), and local angular induction factor (alpha_prime)
   
    lambda_r = (rot_velocity.*radius)/wind_vel;
    lambda = max(lambda_r);
    
    alpha_prime_r = 0.5 * (-1 + sqrt(1 + (4./(lambda_r.^2))*alpha*(1 - alpha)));
end

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
    

    CL_values = zeros(sections, 1); % Empty array to fill with drag coefficeint vlaues
    CD_values = zeros(sections, 1); % Empty array to fill with lift coefficeint vlaues

    % Loop of each section along the blade to select airfoil data and
    % determine coefficients of lift and drag

    for i = 1:sections
        section_airfoil = char(names(i));
        c_AoA = AoA_deg(i); % AoA that is being looped through
    
        % Hard Code CD and CL for cylindrical blade sections using
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
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear', 'extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear', 'extrap');
            
        elseif contains(section_airfoil, 'DU 93-W-210')
            airfoil_table = airfoils.DU93;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear', 'extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear', 'extrap');
            
        elseif contains(section_airfoil, 'DU 96-W-180')
            airfoil_table = airfoils.DU96;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear', 'extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear', 'extrap');
            
        elseif contains(section_airfoil, 'DU 97-W-300')
            airfoil_table = airfoils.DU97;
            CL_values(i) = interp1(airfoil_table.AoA, airfoil_table.CL, c_AoA, 'linear', 'extrap');
            CD_values(i) = interp1(airfoil_table.AoA, airfoil_table.CD, c_AoA, 'linear', 'extrap');
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
    if Re < 2*10^5
        C_D = 11 * Re.^(-0.75) + 0.9 * (1.0 - exp(-1000./Re))...
        + 1.2 * (1.0 -exp(-(Re./4500).^0.7));
    elseif Re <= 5*10^5
        C_D = 10.^(0.32*tanh(44.4504 - 8 * log10 (Re))-0.238793158);
    else
        C_D = 0.1 * log10(Re) - 0.2533429;
    end
end