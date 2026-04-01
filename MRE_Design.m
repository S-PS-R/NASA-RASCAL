%% Design of MSE Reactor (Step 2 of Refinement)
clc
clear
%We need 2 reactors - analysis for 1 reactor using Highlands regolith

T_op = 2100; %electroklysis temperature
M_m = 42.5; %molten regolith mass
gamma = 0.35; % kg O2 per kg regolith
mol_mass_O2 = 0.031999; %kg/mol
Faraday_const = 96485.3321; %C/mol
n_I = 0.82; %current efficiency
T_melt = 1500;

%Volume regolith
p_solid = 1500 %kg/m^3
V_regolith_solid = M_m/p_solid
p_molten = 6.345*10^4/(24.11+0.001206*(T_op-1873))
V_regolith_molten = M_m/p_molten

%Molar mass anorthosite
mol_mass_anorthite = 278.2*10^-3; %kg/mol
n_anorthite = M_m/mol_mass_anorthite;

%%Calculation for average current
C = (M_m*gamma)/mol_mass_O2 * 4*Faraday_const/n_I;
t_batch = 22*3600; %seconds
I = C/(t_batch - 2*60) %A current


%-----Calculation for thermal conductibity of the wall----
%this ensures the reactor can support joule heated cold wall operation
%no molten regolith touches the reactor wall
phi = 1.1;
beta = I^0.215*11*log((T_op - 1525)/(I^0.0105*5387)) - 11*0.351*-5.09 + phi*11*0.343*0.351 + I^0.16*phi*5.31*0.351*log(M_m/(4.99*I^0.435)) - 0.15;
k_wall = (I^0.16*phi*0.351*log((4.99*I^0.435)/M_m) - I^0.215*11*-0.0289*log((T_op - 1525)/(I^0.0105*5387)))/beta

I=I/10^3;

%Diameter Calculations
D = -5.09 - I^0.215/0.351 * log((T_op - 1525)/(I^0.0105*5387))*(1 - 0.0289/(k_wall+0.15))


%Electrode seperation distance
I_avg = I*10^3;
delta_e = 0.0011*I_avg^0.4348*k_wall^0.0570*exp(2.5738*(D + 0.1945)/(I_avg^0.4646*(0.6487/k_wall + 1)^0.3522))*100 %cm

height = 1.5*delta_e*10^-2
V_reactor = pi*D*height

%% Thermoelectric Data
syms T

% Piecewise Specfic Heat Model for Highlands Regolith
%Cp Formula kJ/kg-K (T<350K)
ca = -2.32*10^-2;
cb = 2.13*10^-3;
cc = 1.5*10^-5;
cd = -7.37*10^-8;
ce = 9.66*10^-11;
Cp_350 = (ca + cb*T + cc*T^2 + cd*T^3 + ce*T^4);

%Cp Formula kJ/kg-K (350K < T < 1500)
sum_xi_cf_i = 9.530*10^2;
sum_xi_cg_i = 2.524*10^-1;
sum_xi_ch_i = -2.645*10^7;
Cp_1500 = (sum_xi_cf_i + T*sum_xi_cg_i + T^-2*sum_xi_ch_i)/10^3;

%Cp Formula kJ/kg-K (Molten: T>1500)
sum_xi_cl_i = (1.565*10^3)/10^3;
Cp_molten = sum_xi_cl_i;

mol_mass_CaCl2 = 0.11098; %kg/mol

% Heat capacity for solid CaCl2 (298K to 1045K)
Cp_CaCl2_solid = 3*10^-8*T^3 - 6*10^-5*T^2 + 0.0491*T + 62.55 / (mol_mass_CaCl2 * 10^3);

% Heat capacity for liquid CaCl2 (1045K to 3000K) kJ/kg-k
Cp_CaCl2_liquid = 102.5331/(mol_mass_CaCl2 * 10^3); 

%Heat of formations (kJ/kg)
H_liquid_CaCl2 = -774.09/mol_mass_CaCl2;
H_solid_CaCl2 = -795.8/mol_mass_CaCl2;

%Heat of fusions kJ/kg
H_fusion = 478.6; %Highlands regolith
H_fusion_CaCl2 = H_liquid_CaCl2 - H_solid_CaCl2;

% Piecewise Function for thermal conducitibty (continous doesnt work)

%Thermal conducitivity solid (T < 1500) W/(m-K)
k_solid = (0.001561 + 5.426*10^-11*T^3);

%Thermal conducitivity molten (T > 1500) W/(m-K)
k_molten = (exp(-9.332 + 1.409*10^4/T));

%% Power kJ/s 
%Ptotal = Q_heat + P(gibbs) + Qendo + Q_heat_loss
T_moon = 110;

%Step for heating highlands regolith to molten temp
q1 = double(int(Cp_350, T, T_moon, 350));
q2 = double(int(Cp_1500, T, 350, T_melt));

%Step for heat of fusion of regolith
q3 = H_fusion;

%Step for Heating of CaCl2 to molten state
T_melt_CaCl2 = 1045;
q4 = double(int(Cp_CaCl2_solid, T, T_moon, 1045));

%Step for heat of fusion of CaCl2 
q5 = H_fusion_CaCl2;

%Step for heating regolith to operating temp
q6 = double(int(Cp_molten, T, T_melt, T_op));

%Step for heating CaCl2 to operating temp
q7 = double(int(Cp_CaCl2_liquid, T, T_melt_CaCl2, T_op));

%Total Heat required to get to operating temp
%(NO LONGER USING SALT DUE TO HIGHER POWER FOR HEATING)
Q_heat = M_m*(q1+q2+q3+q6); %+ M_CaCl2*(q4+q5+q7)
kW_heat = Q_heat/(2.25*3600) %assuming 3 hours for full heating


%Assuming highlands regolith is 100% anorthite kJ/mol
Hf_anorthosite = -65.84;
[G_SiO2, Hf_SiO2] = computeDeltaG_H(T_op, -974.5, -0.0453, 0.3612, -974.3, 0.0196);
[G_Al2O3, Hf_Al2O3] = computeDeltaG_H(T_op, -1154, -0.1110, 0.6038, -1154, 0.0481);
[G_CaO, Hf_CaO] = computeDeltaG_H(T_op, -1494, -0.0903, 0.6732, -1491, 0.0377);
[G_MgO, Hf_MgO] = computeDeltaG_H(T_op, -1412, -0.1128, 0.7845, -1409, 0.0479);
[G_FeO, Hf_FeO] = computeDeltaG_H(T_op, -522.9, -0.0130, 0.1561, -522.9, 0.0056);
[G_TiO2, Hf_TiO2] = computeDeltaG_H(T_op, -913.6, -0.0318, 0.2711, -913.5, 0.0137);
[G_Na2O, Hf_Na2O] = computeDeltaG_H(T_op, -1196, -0.2028, 1.218, -1196, 0.0878);
[G_P2O5, Hf_P2O5] = computeDeltaG_H(T_op, -1287, -0.0521, 0.5866, -1286, 0.0221);
[G_K2O, Hf_K2O] = computeDeltaG_H(T_op, -1182, -0.3444, 1.753, -1184, 0.1517);

function [deltaG, deltaH] = computeDeltaG_H(T, a_G, b_G, c_G, a_H, b_H)
    % Compute Delta G(T)
    deltaG = a_G + b_G * T * log10(T) + c_G * T;
    
    % Compute Delta H(T)
    deltaH = a_H + b_H * T;
end


% Constants
R = 8.314/1000; % kJ/(mol·K)
T_op = 2100;
rho_regolith = 6.345*10^4/(24.11+0.001206*(T_op-1873)); %kg/m^3
porosity = 0.83;

% Oxide data: {Name, G (kJ/mol oxide), Molar mass (kg/mol oxide), Mass fraction (kg oxide/kg regolith), Electrons per molecule, Hf (kJ/mol oxide}
oxides = {
    'SiO2',  G_SiO2, 60.08/1000, 0.475, 4, Hf_SiO2;
    'TiO2',  G_TiO2, 79.87/1000, 0.015, 4, Hf_TiO2;
    'CaO',   G_CaO, 56.08/1000, 0.105, 2, Hf_CaO;
    'MgO',   G_MgO, 40.30/1000, 0.09, 2, Hf_MgO;
    'FeO',   G_FeO, 71.84/1000, 0.035, 2, Hf_FeO;
    'K2O',   G_K2O, 94.20/1000, 0.008, 2, Hf_K2O;
    'Na2O',  G_Na2O, 61.98/1000, 0.028, 2, Hf_Na2O;
    'Al2O3', G_Al2O3, 101.96/1000, 0.15, 6, Hf_Al2O3;
    'P2O5',  G_P2O5, 141.94/1000, 0.007, 10, Hf_P2O5;
};

%Step for breaking regolith into oxides kJ
Q8 = 0;
for i = 1:size(oxides, 1)
    name = oxides{i, 1}; % Oxide name
    Hf_oxide = oxides{i,6}; % Heat of formation for oxide in kJ/mol oxide
    m_kg_per_kg_regolith = oxides{i, 4}; % Initial mass of oxide (kg per kg regolith)
    M_kg_per_mol = oxides{i, 3}; % Molar mass of oxide (kg/mol)
    
    % Convert initial mass to moles (moles of oxide per kg regolith)
    moles_per_kg_regolith = m_kg_per_kg_regolith / M_kg_per_mol;
    
    % Calculate the heat absorbed for the oxide: kJ
    H_R = M_m * moles_per_kg_regolith * -Hf_oxide;
    
    % Accumulate total heat absorbed
    Q8 = Q8 + H_R;
end
Q8 = Q8 + n_anorthite*Hf_anorthosite

% Assuming 10 hours for regolith breakdown into oxides (for kW calculation)
kW_oxide = Q8 / (11.75 * 3600) % Conversion from kJ to kW over 10 hours


% Electrochemical Approach using Farday's Law

kW_electrolysis_2 = 0;
V_req_list2 = zeros(size(oxides, 1), 1);
time = 8*3600; %11

% Loop through each oxide to calculate power and minimum voltage
for i = 1:size(oxides, 1)
    G_oxide = oxides{i, 2};
    mol_mass = oxides{i, 3}; % kg/mol oxide
    mass_fraction = oxides{i, 4}; % kg oxide/kg regolith
    n_e = oxides{i, 5}; % Number of electrons per reaction

    % Moles of oxide per kg regolith
    mol_oxide = mass_fraction / mol_mass;

    % Minimum voltage calculation
    V_req = -G_oxide*10^3 / (n_e * Faraday_const);
    V_req_list2(i) = V_req;

    %Current
    current = n_e*mol_oxide*Faraday_const/time;

    %Power
    P_oxide = V_req*current/n_I;
    kW_electrolysis_2 = kW_electrolysis_2 + P_oxide;


    % Display power, current, and voltage for each oxide
    fprintf('%s: Minimum Voltage = %.2f V\n', ...
        oxides{i,1}, V_req);

    
end


avg_G_oxide = (oxides{1, 2}+oxides{2,2}+oxides{4,2}+oxides{5,2}+oxides{8,2})/5;

%Steady state heat sink 
avg_Hf_oxide = (oxides{1, 6}+oxides{2,6}+oxides{4,6}+oxides{5,6}+oxides{8,6})/5;
P_endo = ((M_m*0.35/mol_mass_O2)*(-avg_Hf_oxide+avg_G_oxide))/time
P_electrolysis = ((M_m*0.35/mol_mass_O2)*-avg_G_oxide)/time
P_total_electrolysis = P_endo + P_electrolysis

% Heat loss Q_dot
P_heat_loss = (1.59*10^6*(exp(-8.29/(D+1.08))/(1.15+ 1/(k_wall + 0.0357))) + 222*k_wall)/10^3 %kW

% Total Power
%Ptotal = P_heat + P(gibbs) + Pendo + P_heat_loss
P_total = (kW_heat + kW_oxide) + P_electrolysis + P_endo + P_heat_loss


%% I2R Loss

r_cathode = (0.225)/2;
r_anode = (0.183)/2;

Cathode_SA = 2*pi*(r_cathode)*(0.01+r_cathode);
Anode_SA = 2*pi*(r_anode)*(0.01+r_anode);

Molybendum_epsilon = 0.5;
Iridium_epsilon = 0.3;
sigma = 5.67*10^-8;

P = 15000;

T_anode  = (P/(sigma*Anode_SA*Molybendum_epsilon))^(1/4)
T_cathode = (P/(sigma*Cathode_SA*Iridium_epsilon))^(1/4)




