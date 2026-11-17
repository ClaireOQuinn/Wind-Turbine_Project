%% ========================================================================
%  ME 4053 - Project 2: Wind Turbine Analysis (Clipper Liberty C96)
%  Deliverable 5 (Standalone Script): Structural Analysis of Tower
%  Includes:
%    Part 1: Static Deflection (ODE45) & Bending Stress
%    Part 2: 1-DOF Tower Vibration (pulseFFT) & Tip Superposition
%    Part 3: Global Column Buckling (Euler / Johnson)
%    Part 4: High-Cycle Fatigue Life & Log-Scale S-N Plot
%  Group 29 Operating Point: U_hub = 16.4 m/s, Omega = 15.5 RPM
%% ========================================================================
clc; clear; close all;

%% 1. CONSTANTS & INPUTS FROM DELIVERABLE 4 (GROUP 29)
% Turbine & Operating Specifications
R                 = 48;             % [m] blade radius
hub_height        = 80.4;           % [m] hub height
rotor_speed_range = [9.6, 15.5];    % [RPM] min and max rotor speed
n_blades          = 3;              % [-] number of blades
U_hub_d4          = 16.4;           % [m/s] Group 29 wind speed (Q4 & Q5)
omega_d4          = rotor_speed_range(2) * (pi / 30); % [rad/s] 15.5 RPM
alpha_abl         = 1/7;            % [-] atmospheric boundary layer exponent

% Group 29 Rotor Thrust from Deliverable 4 (at U = 16.4 m/s, pitch = 7.2 deg)
% (When merging into Master Script later, link directly to F_thrust_d4)
F_thrust          = 257.26e3;       % [N] total axial rotor thrust at hub

% Fluid & Material Properties (ASTM A572 Grade 50 Steel)
dens_air             = 1.1;         % [kg/m^3] air density
visc_air             = 1.8e-5;      % [N*s/m^2] dynamic viscosity of air
dens_steel           = 7850;        % [kg/m^3] steel density
strength_tens_steel  = 450e6;       % [Pa] ultimate tensile strength (Su)
strength_yield_steel = 345e6;       % [Pa] yield strength (Sy)
E_steel              = 200e9;       % [Pa] Young's modulus
g                    = 9.81;        % [m/s^2] gravitational acceleration

% Component Masses (Note: m_rotor = 55.7e3 kg already includes hub + 3 blades)
m_rotor   = 55.7e3;                 % [kg] rotor assembly mass
m_nacelle = 94.3e3;                 % [kg] nacelle mass
m_top     = m_rotor + m_nacelle;    % [kg] total Rotor-Nacelle Assembly (RNA) mass

% Import Tower Data
if isfile('towerSpecs.xlsx')
    tower_data = readtable('towerSpecs.xlsx');
else
    tower_data = readtable('towerSpecs.csv');
end

%% 2. PART 1: STATIC DEFORMATION (ODE45) & STATIC NORMAL STRESSES
[z_tower, w_static, w_hub_static, sigma_axial, sigma_bending, sigma_static_max, ...
 SF_yield_static, M_bend, V_shear, P_axial, A_tower, I_tower, OD_tower, m_tower_total] = ...
    analyzeTowerStatic(tower_data, hub_height, U_hub_d4, F_thrust, m_top, ...
    dens_air, visc_air, dens_steel, E_steel, strength_yield_steel, alpha_abl, g);

[peak_static_stress, idx_peak] = max(sigma_static_max);
fprintf('================ DELIVERABLE 5: STRUCTURAL RESULTS ================\n');
fprintf('--- Part 1: Static Deflection & Stress (ODE45) ---\n');
fprintf('Total Tower Mass        : %.1f metric tons\n', m_tower_total / 1000);
fprintf('Base Bending Moment     : %.2f MN*m\n', M_bend(1) / 1e6);
fprintf('Max Static Normal Stress: %.2f MPa (at z = %.2f m)\n', peak_static_stress / 1e6, z_tower(idx_peak));
fprintf('Static Yield Safety Fac.: %.2f (Yield = 345 MPa)\n', SF_yield_static);
fprintf('Tower Top Deflection    : %.1f mm (at z = 77.7 m)\n', w_static(end) * 1000);
fprintf('Hub Center Deflection   : %.1f mm (at z = 80.4 m)\n', w_hub_static * 1000);

%% 3. PART 2: 1-DOF TOWER VIBRATION (pulseFFT) & STRESS SUPERPOSITION
zeta_tower = 0.015;                 % [-] structural damping ratio (1.5% for steel)
nHarmonics = 50;                    % [-] number of Fourier harmonics (< 500)

[t_vib, x_dyn, F_shadow_pulse, k_tower, m_eq, omega_n, f_n, ...
 x_dyn_amp, w_total_max, sigma_total_max, SF_yield_total] = ...
    analyzeTowerVibration(z_tower, OD_tower, I_tower, w_static, sigma_axial, ...
    sigma_bending, R, omega_d4, n_blades, F_thrust, m_top, m_tower_total, ...
    E_steel, strength_yield_steel, zeta_tower, nHarmonics);

fprintf('\n--- Part 2: 1-DOF Vibration & Tip Superposition ---\n');
fprintf('Equivalent Stiffness k  : %.3f MN/m\n', k_tower / 1e6);
fprintf('Natural Frequency wn    : %.3f rad/s (%.3f Hz)\n', omega_n, f_n);
fprintf('3P Blade-Pass Freq.     : %.3f Hz\n', n_blades * rotor_speed_range(2) / 60);
fprintf('Oscillatory Amplitude   : %.2f mm\n', x_dyn_amp * 1000);
fprintf('Peak Dynamic Deflection : %.2f mm (including shadow DC offset)\n', max(abs(x_dyn)) * 1000);
fprintf('Total Tip Deflection    : %.2f mm (Static + Dynamic)\n', w_total_max * 1000);
fprintf('Total Max Bending Stress: %.2f MPa (Superposition at base)\n', sigma_total_max / 1e6);
fprintf('Total Yield Safety Fac. : %.2f (Yield = 345 MPa)\n', SF_yield_total);

%% 4. PART 3: COLUMN BUCKLING ANALYSIS (EULER / JOHNSON)
[L_e, Sr_crit, slend_avg, regime_avg, S_cr_avg, P_cr_avg, SF_buck_avg, ...
 slend_min, regime_min, S_cr_min, P_cr_min, SF_buck_min] = ...
    analyzeTowerBuckling(z_tower, A_tower, I_tower, P_axial, E_steel, strength_yield_steel);

fprintf('\n--- Part 3: Column Buckling Analysis (L_e = 2.1*L) ---\n');
fprintf('Effective Length L_e    : %.2f m\n', L_e);
fprintf('Euler/Johnson Threshold : %.2f\n', Sr_crit);
fprintf('Max Compressive Load    : %.2f MN (at base z = 0 m)\n', P_axial(1) / 1e6);
fprintf('Avg Section (%s)     : L_e/rho = %.2f | S_cr = %.2f MPa | P_cr = %.2f MN | SF = %.2f\n', ...
    regime_avg, slend_avg, S_cr_avg / 1e6, P_cr_avg / 1e6, SF_buck_avg);
fprintf('Top Section (%s)     : L_e/rho = %.2f | S_cr = %.2f MPa | P_cr = %.2f MN | SF = %.2f\n', ...
    regime_min, slend_min, S_cr_min / 1e6, P_cr_min / 1e6, SF_buck_min);

%% 5. PART 4: FATIGUE LIFE ANALYSIS (BENDING ONLY)
% ---> ENTER TEXTBOOK LOOKUP VALUES FOR Kt AND q BELOW <---
% 1) Tower Foundation (Bottom, z = 0 m):
%    Geometric ratios: D/d = 4420/4159 = 1.063,  r/d = 5/4159 = 0.00120
Kt_bot = 2;   % <-- REPLACE NaN with theoretical stress concentration factor Kt (foundation)
q_bot  = 0.67;   % <-- REPLACE NaN with notch sensitivity q (foundation)

% 2) Tower-Nacelle Joint (Top, z = 77.7 m):
%    Geometric ratios: D/d = 3036/3030 = 1.002,  r/d = 5/3030 = 0.00165
Kt_top = 2;   % <-- REPLACE NaN with theoretical stress concentration factor Kt (nacelle)
q_top  = 0.67;   % <-- REPLACE NaN with notch sensitivity q (nacelle)

fatigue_res = analyzeTowerFatigue(z_tower, hub_height, OD_tower, I_tower, M_bend, ...
    R, omega_d4, n_blades, F_thrust, omega_n, zeta_tower, nHarmonics, ...
    strength_tens_steel, Kt_bot, q_bot, Kt_top, q_top);

fprintf('\n--- Part 4: High-Cycle Fatigue Analysis (Bending Only) ---\n');
fprintf('Modified Endurance Limit Se : %.2f MPa\n', fatigue_res.Se_bot / 1e6);
fprintf('Tip Forces (F_mean / F_alt) : %.2f kN / %.2f kN (F_DC = %.2f kN)\n', ...
    fatigue_res.F_mean / 1e3, fatigue_res.F_alt / 1e3, fatigue_res.F_DC / 1e3);
fprintf('Bottom Foundation (z = 0 m) :\n');
fprintf('  Nominal Stress (Mean/Bend): %.2f MPa / %.3f MPa\n', ...
    fatigue_res.sig_mean_bot / 1e6, fatigue_res.sig_bend_bot / 1e6);
fprintf('  Kf = %.3f | sig_alt = %.3f MPa | SF (Se/sig_alt) = %.2f\n', ...
    fatigue_res.Kf_bot, fatigue_res.sig_alt_bot / 1e6, fatigue_res.SF_bot);
fprintf('Top Tower-Nacelle (z = 77.7 m):\n');
fprintf('  Nominal Stress (Mean/Bend): %.2f MPa / %.3f MPa\n', ...
    fatigue_res.sig_mean_top / 1e6, fatigue_res.sig_bend_top / 1e6);
fprintf('  Kf = %.3f | sig_alt = %.3f MPa | SF (Se/sig_alt) = %.2f\n', ...
    fatigue_res.Kf_top, fatigue_res.sig_alt_top / 1e6, fatigue_res.SF_top);
if isnan(Kt_bot) || isnan(q_bot) || isnan(Kt_top) || isnan(q_top)
    fprintf('  [NOTE]: Enter textbook Kt and q values in Section 5 to populate Kf, sig_alt, and SF.\n');
end
fprintf('===================================================================\n');

%% 6. REPORT PLOTS (Q5)
% Plot 1: Static Deflection Plot (Q5)
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

% Plot 3: Dynamic Deflection Plot (Q5)
figure('Name', 'Q5: Dynamic Deflection Plot', 'Color', 'w')
plot(t_vib, x_dyn * 1000, 'b-', 'LineWidth', 1.6)
xlabel('Time t (s)')
ylabel('Dynamic Tip Deflection x_{dyn} (mm)')
title('1-DOF Tower Dynamic Deflection Under Blade Shadow Effect')
grid on

% Plot 4: S-N Plot in Log Scale (Q5 - Required by Rubric)
figure('Name', 'Q5: S-N Plot in Log Scale', 'Color', 'w')
loglog(fatigue_res.N_curve, fatigue_res.Sf_curve / 1e6, 'k-', 'LineWidth', 2.0)
hold on
legend_entries = {'S-N Fatigue Curve (ASTM A572 Gr. 50)'};

if isfinite(fatigue_res.sig_alt_bot)
    yline(fatigue_res.sig_alt_bot / 1e6, 'r--', 'LineWidth', 1.5)
    legend_entries{end+1} = sprintf('Foundation \\sigma_{alt} (%.2f MPa)', fatigue_res.sig_alt_bot / 1e6);
end
if isfinite(fatigue_res.sig_alt_top)
    yline(fatigue_res.sig_alt_top / 1e6, 'b-.', 'LineWidth', 1.5)
    legend_entries{end+1} = sprintf('Nacelle Joint \\sigma_{alt} (%.2f MPa)', fatigue_res.sig_alt_top / 1e6);
end

xlim([1e3, 1e8])
ylim([0.1, 1000])
xlabel('Number of Cycles to Failure, N (cycles)')
ylabel('Fatigue Strength S_f / Alternating Stress \sigma_{alt} (MPa)')
title('Log-Scale S-N Diagram for Wind Turbine Tower')
legend(legend_entries, 'Location', 'southwest')
grid on

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function [C_D] = cylinderCD(Re)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: cylinderCD
%
% PURPOSE:
% Calculates the drag coefficient for a circular cylinder in cross-flow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Re <= 0
        C_D = 0;
    elseif Re < 2e5
        C_D = 11 * Re.^(-0.75) + 0.9 * (1.0 - exp(-1000./Re)) ...
              + 1.2 * (1.0 - exp(-(Re./4500).^0.7));
    elseif Re <= 5e5
        C_D = 10.^(0.32 * tanh(44.4504 - 8 * log10(Re)) - 0.238793158);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_raw    = tower_data{:, 1} / 1000;
    OD_raw   = tower_data{:, 2} / 1000;
    t_raw    = tower_data{:, 3} / 1000;

    L_tower  = max(z_raw);              % [m] steel tower height (77.7 m)
    z_offset = hub_height - L_tower;    % [m] rigid nacelle offset (2.7 m)

    N_nodes = 1000;
    z  = linspace(0, L_tower, N_nodes)';
    OD = interp1(z_raw, OD_raw, z, 'linear');
    t  = interp1(z_raw, t_raw,  z, 'linear');
    ID = OD - 2 .* t;

    A = (pi / 4)  .* (OD.^2 - ID.^2);
    I = (pi / 64) .* (OD.^4 - ID.^4);

    mu_tower      = dens_steel .* A;
    m_tower_total = trapz(z, mu_tower);

    % Atmospheric Boundary Layer wind profile & distributed tower drag
    U_tower  = U_hub .* (z ./ hub_height).^alpha_abl;
    Re_tower = (dens_air .* U_tower .* OD) ./ visc_air;
    CD_tower = arrayfun(@cylinderCD, Re_tower);
    q_wind   = 0.5 * dens_air .* (U_tower.^2) .* CD_tower .* OD;

    % Integrate internal shear V(z), moment M(z), and axial force P(z)
    V_shear = F_thrust + flipud(cumtrapz(flipud(z), flipud(q_wind))) * (-1);
    M_bend  = (F_thrust * z_offset) + flipud(cumtrapz(flipud(z), flipud(V_shear))) * (-1);

    w_tower = mu_tower .* g;
    P_axial = (m_top * g) + flipud(cumtrapz(flipud(z), flipud(w_tower))) * (-1);

    % Solve Euler-Bernoulli beam deflection ODE via ODE45
    beamODE = @(z_eval, Y) [ ...
        Y(2); ...
        interp1(z, M_bend, z_eval, 'linear') / (E_steel * interp1(z, I, z_eval, 'linear')) ...
    ];

    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, Y_sol] = ode45(beamODE, z, [0; 0], opts);

    deflection = Y_sol(:, 1);
    slope      = Y_sol(:, 2);
    w_hub      = deflection(end) + slope(end) * z_offset;

    % Static normal stresses & yield factor of safety
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L_tower = z(end);
    I_avg   = mean(I);
    D_avg   = mean(OD);

    % 1. Equivalent Stiffness, Rayleigh Mass, and Natural Frequency
    k_eq    = (3 * E_steel * I_avg) / (L_tower^3);
    m_eq    = m_top + 0.23 * m_tower;
    omega_n = sqrt(k_eq / m_eq);
    f_n     = omega_n / (2 * pi);

    % 2. Tower Shadow Square-Wave Pulse Parameters
    f_1P        = omega_rotor / (2 * pi);
    freq_3P     = n_blades * f_1P;
    T_3P        = 1 / freq_3P;

    r_eff       = R / 2;
    pulseWidth  = (D_avg / r_eff) / omega_rotor;
    pulseHeight = -F_thrust / n_blades;

    % 3. Fourier Decomposition via pulseFFT
    [fourierMag, fourierPhase] = pulseFFT(freq_3P, pulseWidth, pulseHeight, nHarmonics);

    % 4. Steady-State 1-DOF Harmonic Response
    t       = linspace(0, 4 * T_3P, 2000)';
    x_dc    = fourierMag(1) / k_eq;
    x_osc   = zeros(size(t));
    F_pulse = fourierMag(1) * ones(size(t));

    for idx = 2:nHarmonics
        harm_order = idx - 1;
        w_i        = 2 * pi * harm_order * freq_3P;
        r_i        = w_i / omega_n;

        F_pulse = F_pulse + fourierMag(idx) * cos(w_i * t + fourierPhase(idx));

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

function [L_e, Sr_crit, slend_avg, regime_avg, S_cr_avg, P_cr_avg, SF_buck_avg, ...
          slend_min, regime_min, S_cr_min, P_cr_min, SF_buck_min] = ...
          analyzeTowerBuckling(z, A, I, P_axial, E_steel, sigma_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: analyzeTowerBuckling
%
% PURPOSE:
% Evaluates global column buckling of the wind turbine tower using
% effective length L_e = 2.1*L, radius of gyration rho = sqrt(I/A),
% and Euler / Johnson critical buckling stress criteria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L_tower = z(end);                   % [m] structural steel tower height (77.7 m)
    L_e     = 2.1 * L_tower;            % [m] effective column length (fixed-free)
    P_max   = P_axial(1);               % [N] maximum compressive force at base

    % Euler-Johnson transition slenderness ratio
    Sr_crit = sqrt((2 * pi^2 * E_steel) / sigma_y);

    % 1. Average cross-section evaluation
    A_avg     = mean(A);
    I_avg     = mean(I);
    rho_avg   = sqrt(I_avg / A_avg);
    slend_avg = L_e / rho_avg;

    if slend_avg > Sr_crit
        regime_avg = 'Euler';
        S_cr_avg   = (pi^2 * E_steel) / (slend_avg^2);
    else
        regime_avg = 'Johnson';
        S_cr_avg   = sigma_y - ((sigma_y^2) * (slend_avg^2)) / (4 * pi^2 * E_steel);
    end
    P_cr_avg    = S_cr_avg * A_avg;
    SF_buck_avg = P_cr_avg / P_max;

    % 2. Conservative minimum cross-section evaluation (top flange)
    A_min     = A(end);
    I_min     = I(end);
    rho_min   = sqrt(I_min / A_min);
    slend_min = L_e / rho_min;

    if slend_min > Sr_crit
        regime_min = 'Euler';
        S_cr_min   = (pi^2 * E_steel) / (slend_min^2);
    else
        regime_min = 'Johnson';
        S_cr_min   = sigma_y - ((sigma_y^2) * (slend_min^2)) / (4 * pi^2 * E_steel);
    end
    P_cr_min    = S_cr_min * A_min;
    SF_buck_min = P_cr_min / P_max;
end

function res = analyzeTowerFatigue(z, hub_height, OD, I, M_bend, ...
    R, omega_rotor, n_blades, F_thrust, omega_n, zeta, nHarmonics, ...
    Su, Kt_bot, q_bot, Kt_top, q_top)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: analyzeTowerFatigue
%
% PURPOSE:
% Performs bending-only high-cycle fatigue analysis at the tower foundation
% (bottom) and tower-nacelle joint (top) using Norton correction factors,
% pulseFFT harmonic transmissibility (TR), stress concentration (Kf), and
% generates the log-scale S-N curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L_tower  = z(end);
    z_offset = hub_height - L_tower;    % [m] 2.7 m rigid nacelle offset
    D_avg    = mean(OD);

    % 1. Endurance Limit Correction Factors (Bending Only)
    C_load    = 1.0;                    % pure bending (axial & shear neglected)
    C_surface = 0.8;                    % machined surface
    C_temp    = 1.0;                    % ambient temperature
    C_reliab  = 0.897;                  % 90% reliability
    Se_prime  = 0.5 * Su;               % [Pa] uncorrected endurance limit

    % Piecewise C_size based on outer diameter d [mm]
    d_bot_mm = OD(1) * 1000;
    d_top_mm = OD(end) * 1000;
    C_size_bot = calcSizeFactor(d_bot_mm);
    C_size_top = calcSizeFactor(d_top_mm);

    res.Se_bot = C_load * C_size_bot * C_surface * C_temp * C_reliab * Se_prime;
    res.Se_top = C_load * C_size_top * C_surface * C_temp * C_reliab * Se_prime;

    % 2. Dynamic Tip Forces from pulseFFT & Transmissibility Ratio (TR)
    f_1P        = omega_rotor / (2 * pi);
    freq_3P     = n_blades * f_1P;
    T_3P        = 1 / freq_3P;
    r_eff       = R / 2;
    pulseWidth  = (D_avg / r_eff) / omega_rotor;
    pulseHeight = -F_thrust / n_blades;

    [fourierMag, fourierPhase] = pulseFFT(freq_3P, pulseWidth, pulseHeight, nHarmonics);

    res.F_DC   = fourierMag(1);
    res.F_mean = F_thrust + res.F_DC;

    % Reconstruct dynamic tip force F_dyn(t) = F_DC + sum((F_i * TR_i)*cos(...))
    t     = linspace(0, 4 * T_3P, 2000)';
    F_dyn = res.F_DC * ones(size(t));
    for idx = 2:nHarmonics
        harm_order = idx - 1;
        w_i        = 2 * pi * harm_order * freq_3P;
        r_i        = w_i / omega_n;
        TR_i       = 1 / sqrt((1 - r_i^2)^2 + (2 * zeta * r_i)^2);
        psi_i      = atan2(2 * zeta * r_i, 1 - r_i^2);

        F_dyn = F_dyn + (fourierMag(idx) * TR_i) * cos(w_i * t + fourierPhase(idx) - psi_i);
    end

    res.F_alt = (max(F_dyn) - min(F_dyn)) / 2;
    F_max     = res.F_mean + res.F_alt;
    F_min     = res.F_mean - res.F_alt;

    % 3. Bending Moments & Nominal Stresses at Bottom (z = 0) and Top (z = 77.7 m)
    % Retain distributed tower wind-drag moment at the base from Part 1
    M_wind_base = M_bend(1) - F_thrust * hub_height;

    M_max_bot = F_max * hub_height + M_wind_base;
    M_min_bot = F_min * hub_height + M_wind_base;
    sig_max_bot = (M_max_bot * (OD(1) / 2)) / I(1);
    sig_min_bot = (M_min_bot * (OD(1) / 2)) / I(1);
    res.sig_mean_bot = (sig_max_bot + sig_min_bot) / 2;
    res.sig_bend_bot = (sig_max_bot - sig_min_bot) / 2;

    M_max_top = F_max * z_offset;
    M_min_top = F_min * z_offset;
    sig_max_top = (M_max_top * (OD(end) / 2)) / I(end);
    sig_min_top = (M_min_top * (OD(end) / 2)) / I(end);
    res.sig_mean_top = (sig_max_top + sig_min_top) / 2;
    res.sig_bend_top = (sig_max_top - sig_min_top) / 2;

    % 4. Fatigue Stress Concentration Factor Kf & Safety Factor SF
    res.Kf_bot      = 1 + q_bot * (Kt_bot - 1);
    res.sig_alt_bot = res.Kf_bot * res.sig_bend_bot;
    res.SF_bot      = res.Se_bot / res.sig_alt_bot;

    res.Kf_top      = 1 + q_top * (Kt_top - 1);
    res.sig_alt_top = res.Kf_top * res.sig_bend_top;
    res.SF_top      = res.Se_top / res.sig_alt_top;

    % 5. Log-Scale S-N Curve (10^3 <= N <= 10^6, plus endurance plateau to 10^8)
    N_finite   = logspace(3, 6, 200);
    log10_Sf   = (1/3) * log10(res.Se_bot / (0.9 * Su)) .* log10(N_finite) ...
                 + log10(((0.9 * Su)^2) / res.Se_bot);
    Sf_finite  = 10.^log10_Sf;

    N_plateau  = logspace(6, 8, 50);
    Sf_plateau = res.Se_bot * ones(size(N_plateau));

    res.N_curve  = [N_finite, N_plateau(2:end)];
    res.Sf_curve = [Sf_finite, Sf_plateau(2:end)];
end

function C_size = calcSizeFactor(d_mm)
% Piecewise size factor C_size from diameter d in mm
    if d_mm <= 8
        C_size = 1.0;
    elseif d_mm <= 250
        C_size = 1.189 * (d_mm ^ -0.097);
    else
        C_size = 0.6;
    end
end
