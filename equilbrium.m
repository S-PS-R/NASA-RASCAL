%% This section assumes thermodynamic equilbrium which apprently isn't true

% Preallocate arrays for storing concentrations
num_oxides = size(oxides, 1);
C_initial_array = zeros(num_oxides, 1);
C_eq_ion_array = zeros(num_oxides, 1);
C_eq_2_array = zeros(num_oxides, 1);

% Loop through each oxide and calculate equilibrium properties
for i = 1:2
    name = oxides{i,1};
    G_O2_kJ = oxides{i,2}; % Gibbs free energy in kJ/mol O2
    M_kg_per_mol = oxides{i,3}; % Molar mass (kg/mol)
    m_kg_per_kg_regolith = oxides{i,4}; % Initial mass (kg oxide per kg regolith)

    % Convert Gibbs free energy from kJ/mol O2 to kJ/mol of oxide
    G_oxide_J = G_O2_kJ * stoichiometry.(name);

    % Calculate the equilibrium constant
    K = exp(-G_oxide_J / (R * T_op));

    % Convert initial mass to molar concentration [mol/m^3]
    a = (m_kg_per_kg_regolith / M_kg_per_mol) * rho_regolith;

    % Assume reaction: Oxide ⇌ Ions
    syms x % Concentration of ions
    c = a - x; % Remaining concentration of Oxide
    eqn = K == (x*(2*x)^2)/c;  % Equilibrium equation
    sol = double(vpasolve(eqn, x)); % Solve for x
    
    % Select the valid solution (positive and real)
    sol = sol(sol > 0); % Keep only positive solutions

    % Store results in arrays
    C_initial_array(i) = a;
    C_eq_ion_array(i) = sol;

     % Display results
    fprintf('%s: K = %.2e, C_initial = %.4f M, C_eq_ion = %.4f M\n', ...
            name, K, a, sol(1));
end

for i = 3:5
    name = oxides{i,1};
    G_O2_kJ = oxides{i,2}; % Gibbs free energy in kJ/mol O2
    M_kg_per_mol = oxides{i,3}; % Molar mass (kg/mol)
    m_kg_per_kg_regolith = oxides{i,4}; % Initial mass (kg oxide per kg regolith)

    % Convert Gibbs free energy from kJ/mol O2 to kJ/mol of oxide
    G_oxide_J = G_O2_kJ * stoichiometry.(name);

    % Calculate the equilibrium constant
    K = exp(-G_oxide_J / (R * T_op));

    % Convert initial mass to molar concentration [mol/m^3]
    a = (m_kg_per_kg_regolith / M_kg_per_mol) * (rho_regolith / (1 - porosity));

    % Assume reaction: Oxide ⇌ Ions
    syms x % Concentration of ions
    c = a - x; % Remaining concentration of Oxide
    eqn = K == x^2/c;  % Equilibrium equation
    sol = double(vpasolve(eqn, x)); % Solve for x
    
    % Select the valid solution (positive and real)
    sol = sol(sol > 0); % Keep only positive solutions

    % Store results in arrays
    C_initial_array(i) = a;
    C_eq_ion_array(i) = sol;

     % Display results
    fprintf('%s: K = %.2e, C_initial = %.4f M, C_eq_ion = %.4f M\n', ...
            name, K, a, sol(1));
end

for i = 6:7
    name = oxides{i,1};
    G_O2_kJ = oxides{i,2}; % Gibbs free energy in kJ/mol O2
    M_kg_per_mol = oxides{i,3}; % Molar mass (kg/mol)
    m_kg_per_kg_regolith = oxides{i,4}; % Initial mass (kg oxide per kg regolith)

    % Convert Gibbs free energy from kJ/mol O2 to kJ/mol of oxide
    G_oxide_J = G_O2_kJ * stoichiometry.(name);

    % Calculate the equilibrium constant
    K = exp(-G_oxide_J / (R * T_op));

    % Convert initial mass to molar concentration [mol/m^3]
    a = (m_kg_per_kg_regolith / M_kg_per_mol) * (rho_regolith / (1 - porosity));

    % Assume reaction: Oxide ⇌ Ions
    syms x % Concentration of ions
    c = a - x; % Remaining concentration of Oxide
    eqn = K == 4*x^3/c;  % Equilibrium equation
    sol = double(vpasolve(eqn, x)); % Solve for x
    
    % Select the valid solution (positive and real)
    sol = sol(sol > 0); % Keep only positive solutions

    % Store results in arrays
    C_initial_array(i) = a;
    C_eq_ion_array(i) = sol;

     % Display results
    fprintf('%s: K = %.2e, C_initial = %.4f M, C_eq_ion = %.4f M\n', ...
            name, K, a, sol(1));
end

for i = 8:8
    name = oxides{i,1};
    G_O2_kJ = oxides{i,2}; % Gibbs free energy in kJ/mol O2
    M_kg_per_mol = oxides{i,3}; % Molar mass (kg/mol)
    m_kg_per_kg_regolith = oxides{i,4}; % Initial mass (kg oxide per kg regolith)

    % Convert Gibbs free energy from kJ/mol O2 to kJ/mol of oxide
    G_oxide_J = G_O2_kJ * stoichiometry.(name);

    % Calculate the equilibrium constant
    K = exp(-G_oxide_J / (R * T_op));

    % Convert initial mass to molar concentration [mol/m^3]
    a = (m_kg_per_kg_regolith / M_kg_per_mol) * (rho_regolith / (1 - porosity));

    % Assume reaction: Oxide ⇌ Ions
    syms x % Concentration of ions
    c = a - x; % Remaining concentration of Oxide
    eqn = K == 108*x^5/c;  % Equilibrium equation
    sol = double(vpasolve(eqn, x)); % Solve for x
    
    % Select the valid solution (positive and real)
    sol = sol(sol > 0); % Keep only positive solutions

    % Store results in arrays
    C_initial_array(i) = a;
    C_eq_ion_array(i) = sol(1);

    % Display results
    fprintf('%s: K = %.2e, C_initial = %.4f M, C_eq_ion = %.4f M\n', ...
            name, K, a, sol(1));
end

for i = 9:9
    name = oxides{i,1};
    G_O2_kJ = oxides{i,2}; % Gibbs free energy in kJ/mol O2
    M_kg_per_mol = oxides{i,3}; % Molar mass (kg/mol)
    m_kg_per_kg_regolith = oxides{i,4}; % Initial mass (kg oxide per kg regolith)

    % Convert Gibbs free energy from kJ/mol O2 to kJ/mol of oxide
    G_oxide_J = G_O2_kJ * stoichiometry.(name);

    % Calculate the equilibrium constant
    K = exp(-G_oxide_J / (R * T_op));

    % Convert initial mass to molar concentration [mol/m^3]
    a = (m_kg_per_kg_regolith / M_kg_per_mol) * (rho_regolith / (1 - porosity));

    % Assume reaction: Oxide ⇌ Ions
    syms x % Concentration of ions
    c = a - x; % Remaining concentration of Oxide
    eqn = K == 12500*x^7/c;  % Equilibrium equation
    sol = double(vpasolve(eqn, x)); % Solve for x
    
    % Select the valid solution (positive and real)
    sol = sol(sol > 0); % Keep only positive solutions

    % Store results in arrays
    C_initial_array(i) = a;
    C_eq_ion_array(i) = sol(1);

    % Display results
    fprintf('%s: K = %.2e, C_initial = %.4f M, C_eq_ion = %.4f M\n', ...
            name, K, a, sol(1));

end

%Convert concentrations back to mol
% Preallocate arrays for storing equilibrium moles
moles_b = zeros(num_oxides, 1);
moles_c = zeros(num_oxides, 1);
moles_a = zeros(num_oxides, 1);

% Loop through each oxide to calculate equilibrium moles
for i = 1:num_oxides
    % Extract equilibrium concentrations
    a = C_initial_array(i);
    moles_a(i) = a * V_regolith_molten;
    c = C_eq_ion_array(i);

    % Calculate equilibrium moles of oxide left (mol = conc * volume)
    moles_b(i) = (a-c) * V_regolith_molten;

    % Calculate moles of metal/oxygen removed (difference in initial and equilibrium)
    moles_c(i) = c * V_regolith_molten;
end


%Electrolysis energy demand

% Oxides concentration data: name, C_initial, C_eq_ion, deltaG (Gibbs free energy change)
oxides_concentration = {
    'SiO2', 20572.8567, 20570.8244, G_SiO2; % SiO2: C_initial, C_eq_ion, deltaG
    'TiO2', 488.6957, 488.6957, G_TiO2;   % TiO2: C_initial, C_eq_ion, deltaG
    'CaO', 28659.1489, 28657.9448, G_CaO;  % CaO: C_initial, C_eq_ion, deltaG
    'MgO', 34183.7306, 34023.7032, G_MgO;  % MgO: C_initial, C_eq_ion, deltaG
    'FeO', 7457.3361, 3679.6992, G_FeO;    % FeO: C_initial, C_eq_ion, deltaG
    'K2O', 1299.9333, 2.7337, G_K2O;       % K2O: C_initial, C_eq_ion, deltaG
    'Na2O', 6914.9406, 19.8988, G_Na2O;    % Na2O: C_initial, C_eq_ion, deltaG
    'Al2O3', 22518.7058, 22518.5631, G_Al2O3; % Al2O3: C_initial, C_eq_ion, deltaG
    'P2O5', 754.8753, 754.8591, G_P2O5;    % P2O5: C_initial, C_eq_ion, deltaG
};

% Preallocate arrays
moles_ions = zeros(size(oxides_concentration, 1), 1); % Moles of ions produced
energy_demand = zeros(size(oxides_concentration, 1), 1); % Energy demand for each oxide

% Loop through each oxide
for i = 1:size(oxides_concentration, 1)
    
    name = oxides_concentration{i, 1};
    K = oxides_concentration{i, 2}; % Equilibrium constant
    C_initial = oxides_concentration{i, 2}; % Initial concentration of oxide
    C_eq_ion = oxides_concentration{i, 3}; % Equilibrium concentration of ions
    deltaG = oxides_concentration{i, 4}; % Gibbs free energy change for electrolysis

    % Compute Gibbs free energy for the electrolysis reaction (simplified)
    % We assume a simplified form here, as exact calculations will require 
    % the Gibbs free energies of the ions and oxides at the operating temperature.
    % This is a placeholder and would be calculated based on real data.
    % For simplicity, we assume delta G for the electrolysis reaction is proportional to K.
   
    
    % Calculate moles of ions from equilibrium concentration
    moles_ions(i) = C_eq_ion * V_regolith_molten;

    % Calculate electrolysis energy demand (in Joules)
    energy_demand(i) = ((delta_G / Faraday_const) * moles_ions(i))/10^3; % kJ
end

% Total energy demand for all oxides
electrolysis_energy = sum(energy_demand);
kW_electrolysis = electrolysis_energy/(11*3600) %Assumes 11 hours 
