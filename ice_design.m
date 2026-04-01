%Ice regolith
clc
clear

%Blue moon propellant ratio
O2_per_H2 = 5;

%Ice kg per kg regolith (range)
n_ice = 0.05;

%Water extraction efffiency
n_water = 0.55;

%Yearly production demand of H2 (kg/yr)
mass_H2 = 1700;

%Daily production demand of H2 (kg/day)
mass_H2_day = mass_H2/365;
mol_mass_H2 = 0.002; %kg/mol
mol_H2_day = mass_H2_day/mol_mass_H2;

%Moles of H2 = Moles of H20
mol_mass_H20 = 0.0180153; %kg/mol
mass_H20_day = mol_H2_day*mol_mass_H20

%Given extraction effiency, what is icy regolith required (mt/day)
M_regolith_day = ((mass_H20_day/n_ice)/n_water)/1000
mass_regolith_year = M_regolith_day*365


%Water containement
time = 60*60;
m_per_batch = 20; %kg
kg_per_day = (24*3600)*m_per_batch/time
num_rovers = M_regolith_day*10^3/kg_per_day