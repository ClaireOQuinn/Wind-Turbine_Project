%% Wind Turbine Project

clc
clear
clear all

%% Import all Data and Variables

rated_power = 2.5 * (10^6); % [MW] rated turbine power
blade_rad = 48; % [m] blade radius
hub_height = 80.4; % [m] 
cut_in_speed = 4; % [m/s]
rated_speed = 11; % [m/s]
cut_out_speed = 25; % [m/s]
rotor_speed_min = 9.6; % [RPM]
rotor_speed_max = 15.5; % [RPM]
dens_air = 1.1; % [kg/m^3] density of air
visc_air = 8 * (10^-5); % [Ns/m^2] viscosity of air
dens_steel = 7850; % [kg/m^3] density of steel
strength_tens_steel = 450 * (10^6); % [Pa] tensile strength of steel
strength_yield_steel = 345 * (10^6); % [Pa] yield strength of steel

% Import Turbine Data
Blade_data = readtable('BladeProfile.csv'); % Blade data
DU91_data = readtable('DU91-W2-250.csv'); % Airfoil data
DU93_data = readtable('DU93-W-210.csv'); % Airfoil data
DU96_data = readtable('DU96-W-180.csv'); % Airfoil data
DU97_data = readtable('DU97-W-300.csv'); % Airfoil data
tower_data = readtable('towerSpecs.csv'); % Tower data
Blade_data.DistanceFromCenterOfRotation = Blade_data.DistanceFromCenterOfRotation/1000; %converting to m

%% Deliverable 1: Determining coefficient of power and coefficient of thrust
wind_vel = 10; % [m/s]
rot_velocity = 14; % [RPM]
pitch = 0; % [degrees]



