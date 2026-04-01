%% Step 2a: Heating rock to melting points calculation
syms T

Cp_anorthosite = 1*10^-10*T^3 + 4*10^-8*T^2 + 0.0006*T + 0.6217;
Cp_basalt = -1*10^-10*T^3 + 3*10^-7*T^2 + 7*10^-5*T + 0.7436;

%Step 2a Temp range for anorthosite to melt (k)
T_moon = 110; 
T_melt_anorthosite = 1826.15;
T_melt_basalt = 1770;

% Compute the definite integral (kJ/kg)
specific_energy_anorthsite2a = double(int(Cp_anorthosite, T, T_moon, T_melt_anorthosite))
specific_energy_basalt2a = double(int(Cp_basalt, T, T_moon, T_melt_basalt))


%% Step 2b: Heat of fusion calculation

%Change dHf for anorthosite from kcal/mol to kJ/kg
Hf_anorthosite_kcal_mol = 32.4;  %[kcal/mol]
kcal_to_kJ = 4.184;  % [kJ/kcal]
molar_mass_anorthite = 278.16/1000;  % Molar mass of CaAl2Si2O8 in kg/mol

% Convert kcal/mol to kJ/kg
Hf_anorthosite = Hf_anorthosite_kcal_mol * kcal_to_kJ / molar_mass_anorthite


%Hf formula for basalt
Hf_basalt = -118.48+0.36311*T_melt_basalt

%% Step 2c: Heating melted rock to operating temp

%Operating MSE temperature
T_MSE = 1850; %K
specific_energy_anorthsite2c = double(int(Cp_anorthosite, T, T_melt_anorthosite, T_MSE))
specific_energy_basalt2c = double(int(Cp_basalt, T, T_melt_basalt, T_MSE))

%% Step 2d: Heat CaCl2 to melting -> melting -> operating temp

T_solid = linspace(0, 1045, 500); % Solid phase from 0 to 1045 K
T_liquid = linspace(1045, 2500, 500); % Liquid phase from 1045 to 2500 K

% Convert to reduced temperature t = T/1000
t_solid = T_solid / 1000;
t_liquid = T_liquid / 1000;

% Shomate equation coefficients for solid phase (298-1045K)
A_s = 87.29581;
B_s = -35.07644;
C_s = 44.12781;
D_s = -9.848049;
E_s = -0.674264;

% Compute Cp for solid phase
Cp_solid = A_s + B_s.*t_solid + C_s.*t_solid.^2 + D_s.*t_solid.^3 + E_s./t_solid.^2;

% Shomate equation coefficients for liquid phase (1045-3000K)
A_l = 102.5331;
B_l = 2.006254e-8;
C_l = -9.234423e-10;
D_l = 1.407022e-9;
E_l = 2.609206e-9;

% Compute Cp for liquid phase
Cp_liquid = A_l + B_l.*t_liquid + C_l.*t_liquid.^2 + D_l.*t_liquid.^3 + E_l./t_liquid.^2;

% Plot the results
figure;
hold on;
plot(T_solid, Cp_solid, 'b', 'LineWidth', 2); % Solid phase in blue
plot(T_liquid, Cp_liquid, 'r', 'LineWidth', 2); % Liquid phase in red

% Formatting the plot
xlabel('Temperature (K)', 'FontSize', 12);
ylabel('Heat Capacity (J/mol-K)', 'FontSize', 12);
title('Specific Heat of CaCl_2 using Shomate Equation', 'FontSize', 14);
xlim([298, 2500]);
ylim([60, 120]);
legend('Solid Phase', 'Liquid Phase', 'Location', 'NorthWest');
grid on;
hold off;


T_melt_CaCl2 = 1045;
CaCl2_mass = 0.11098; %kg/mol
Cp_CaCl2_298 = A_s + B_s*T/1000 + C_s*(T/1000)^2 + D_s*(T/1000)^3 + E_s/(T/1000)^2;
Cp_CaCl2_1045 = A_l + B_l*T/1000 + C_l*(T/1000)^2 + D_l*(T/1000)^3 + E_l/(T/1000)^2;

M_m = 42.5;
M_CaCl2 = M_m*1

%kJ/kg
specific_energy_CaCl2 = double(int(Cp_CaCl2_298, T, T_moon, T_melt_CaCl2))/(CaCl2_mass*10^3)
%Melting
Hf_CaCl2 = 43.4/CaCl2_mass %not totally accurate since couldn't find information on fusion for CaCl2 only CaCl2*6H20
specific_energy_CaCl2_2 = double(int(Cp_CaCl2_1045, T, T_melt_CaCl2, T_MSE))/(CaCl2_mass *10^3);

P_CaCl2 = (M_CaCl2*(specific_energy_CaCl2+specific_energy_CaCl2_2))/(2.25*3600)

 
%% Step 3: Dissociation into oxides

%Anorthite Calculation all kJ/mol
% [s] = solid, [l] = liquid for if 1850K is low/higher than oxide melting point
Hf_anorthosite = -65.84; %this value was given for 1650K [s]
Hf_CaO = -635.5; % [s]
Hf_Al2O3 = -1669.8; % [s]
Hf_SiO2 = -859.4; % [s]

Hr = Hf_CaO + Hf_Al2O3 + Hf_SiO2 - Hf_anorthosite


