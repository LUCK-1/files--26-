function params = MaterialDatabase(materialChoice)

    switch materialChoice
        case 'LowCarbonSteel'
            params.materialName = 'Low Carbon Steel (AISI 1010)';
            params.b = 2.48e-10;
            params.mu = 80e9;
            params.Qact = 280e3;
            params.R = 8.314;
            params.k1 = 1e9;
            params.k2_0 = 20;
            params.m_k2 = 0.1;
            params.Qk2 = 100e3;
            params.rho0 = 1e10;
            params.rho_critical = 1e13;
            params.rho_saturation = 1e15;
            params.Qnuc = 120e3;
            params.C_nuc = 0.5;
            params.nuc_probability_base = 0.15;
            params.Qgb = 120e3;
            params.Mgb0 = 1e-5;
            params.gamma_gb = 0.5;
            params.gamma_lagb = 0.15;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.3;
            params.maxGrowthFraction = 0.02;
            
        case 'AISI304'
            params.materialName = 'Austenitic Stainless Steel (AISI 304)';
            params.b = 2.58e-10;
            params.mu = 77e9;
            params.Qact = 400e3;
            params.R = 8.314;
            params.k1 = 1.2e9;
            params.k2_0 = 18;
            params.m_k2 = 0.08;
            params.Qk2 = 120e3;
            params.rho0 = 1e10;
            params.rho_critical = 1.5e13;
            params.rho_saturation = 1.2e15;
            params.Qnuc = 150e3;
            params.C_nuc = 0.4;
            params.nuc_probability_base = 0.12;
            params.Qgb = 140e3;
            params.Mgb0 = 8e-6;
            params.gamma_gb = 0.65;
            params.gamma_lagb = 0.2;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.32;
            params.maxGrowthFraction = 0.018;
            
        case 'Aluminum6061'
            params.materialName = 'Aluminum Alloy 6061';
            params.b = 2.86e-10;
            params.mu = 26e9;
            params.Qact = 156e3;
            params.R = 8.314;
            params.k1 = 8e8;
            params.k2_0 = 25;
            params.m_k2 = 0.12;
            params.Qk2 = 80e3;
            params.rho0 = 1e10;
            params.rho_critical = 8e12;
            params.rho_saturation = 5e14;
            params.Qnuc = 100e3;
            params.C_nuc = 0.6;
            params.nuc_probability_base = 0.2;
            params.Qgb = 100e3;
            params.Mgb0 = 2e-5;
            params.gamma_gb = 0.35;
            params.gamma_lagb = 0.10;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.25;
            params.maxGrowthFraction = 0.025;
            
        case 'Copper'
            params.materialName = 'Pure Copper (OFHC)';
            params.b = 2.56e-10;
            params.mu = 48e9;
            params.Qact = 197e3;
            params.R = 8.314;
            params.k1 = 9e8;
            params.k2_0 = 22;
            params.m_k2 = 0.1;
            params.Qk2 = 90e3;
            params.rho0 = 1e10;
            params.rho_critical = 9e12;
            params.rho_saturation = 8e14;
            params.Qnuc = 110e3;
            params.C_nuc = 0.5;
            params.nuc_probability_base = 0.15;
            params.Qgb = 110e3;
            params.Mgb0 = 1.5e-5;
            params.gamma_gb = 0.45;
            params.gamma_lagb = 0.12;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.28;
            params.maxGrowthFraction = 0.02;
            
        case 'Titanium'
            params.materialName = 'Titanium (Ti-6Al-4V)';
            params.b = 2.95e-10;
            params.mu = 44e9;
            params.Qact = 330e3;
            params.R = 8.314;
            params.k1 = 1.1e9;
            params.k2_0 = 18;
            params.m_k2 = 0.08;
            params.Qk2 = 130e3;
            params.rho0 = 1e10;
            params.rho_critical = 1.2e13;
            params.rho_saturation = 1e15;
            params.Qnuc = 140e3;
            params.C_nuc = 0.35;
            params.nuc_probability_base = 0.1;
            params.Qgb = 140e3;
            params.Mgb0 = 6e-6;
            params.gamma_gb = 0.55;
            params.gamma_lagb = 0.16;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.32;
            params.maxGrowthFraction = 0.015;
            
        case 'Nickel718'
            params.materialName = 'Nickel Superalloy (Inconel 718)';
            params.b = 2.54e-10;
            params.mu = 82e9;
            params.Qact = 450e3;
            params.R = 8.314;
            params.k1 = 1.3e9;
            params.k2_0 = 15;
            params.m_k2 = 0.06;
            params.Qk2 = 150e3;
            params.rho0 = 1e10;
            params.rho_critical = 1.8e13;
            params.rho_saturation = 1.5e15;
            params.Qnuc = 170e3;
            params.C_nuc = 0.3;
            params.nuc_probability_base = 0.08;
            params.Qgb = 160e3;
            params.Mgb0 = 5e-6;
            params.gamma_gb = 0.7;
            params.gamma_lagb = 0.20;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.35;
            params.maxGrowthFraction = 0.012;
            
        case 'Magnesium'
            params.materialName = 'Magnesium Alloy (AZ31)';
            params.b = 3.21e-10;
            params.mu = 17e9;
            params.Qact = 135e3;
            params.R = 8.314;
            params.k1 = 7e8;
            params.k2_0 = 28;
            params.m_k2 = 0.15;
            params.Qk2 = 70e3;
            params.rho0 = 1e10;
            params.rho_critical = 6e12;
            params.rho_saturation = 4e14;
            params.Qnuc = 90e3;
            params.C_nuc = 0.7;
            params.nuc_probability_base = 0.25;
            params.Qgb = 90e3;
            params.Mgb0 = 3e-5;
            params.gamma_gb = 0.30;
            params.gamma_lagb = 0.08;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.22;
            params.maxGrowthFraction = 0.03;
            
        case 'Custom'
            params.materialName = 'Custom Material';
            params.b = 2.5e-10;
            params.mu = 50e9;
            params.Qact = 200e3;
            params.R = 8.314;
            params.k1 = 9e8;
            params.k2_0 = 22;
            params.m_k2 = 0.1;
            params.Qk2 = 100e3;
            params.rho0 = 1e10;
            params.rho_critical = 1e13;
            params.rho_saturation = 8e14;
            params.Qnuc = 120e3;
            params.C_nuc = 0.5;
            params.nuc_probability_base = 0.15;
            params.Qgb = 120e3;
            params.Mgb0 = 1e-5;
            params.gamma_gb = 0.5;
            params.gamma_lagb = 0.15;
            params.theta_HAGB = 15;
            params.alpha_taylor = 0.3;
            params.maxGrowthFraction = 0.02;
            
        otherwise
            error('Unknown material: %s\nAvailable: LowCarbonSteel, AISI304, Aluminum6061, Copper, Titanium, Nickel718, Magnesium, Custom', materialChoice);
    end
    
    params.strainIncrement = 0.001;
    params.outputInterval = 20;
    params.gridSize = 300;
    params.cellSize = 1e-6;
    params.initialGrainNum = 40;
    params.initialGrainSize = 50e-6;
    params.neighborhoodType = 'Moore';
    params.maxGrainID = 10000;
    params.deformationType = 'compression';
    params.deformationInterval = 50;
end

function PrintMaterialParameters(params)
    fprintf('\n=============== Material Parameters ===============\n');
    fprintf('Material: %s\n', params.materialName);
    fprintf('---------------------------------------------------\n');
    fprintf('Burgers vector (b): %.2e m\n', params.b);
    fprintf('Shear modulus (mu): %.2e Pa\n', params.mu);
    fprintf('Activation energy (Qact): %.0f kJ/mol\n', params.Qact/1000);
    fprintf('---------------------------------------------------\n');
    fprintf('Dislocation Hardening (k1): %.2e m^-1\n', params.k1);
    fprintf('Dynamic Recovery (k2_0): %.1f\n', params.k2_0);
    fprintf('Strain rate sensitivity (m_k2): %.3f\n', params.m_k2);
    fprintf('Recovery activation energy (Qk2): %.0f kJ/mol\n', params.Qk2/1000);
    fprintf('---------------------------------------------------\n');
    fprintf('Initial dislocation density: %.2e m^-2\n', params.rho0);
    fprintf('Critical dislocation density: %.2e m^-2\n', params.rho_critical);
    fprintf('---------------------------------------------------\n');
    fprintf('Nucleation activation energy: %.0f kJ/mol\n', params.Qnuc/1000);
    fprintf('GB migration activation energy: %.0f kJ/mol\n', params.Qgb/1000);
    fprintf('GB mobility prefactor: %.2e m^4/(J*s)\n', params.Mgb0);
    fprintf('HAGB energy: %.3f J/m^2\n', params.gamma_gb);
    fprintf('LAGB energy coefficient: %.3f J/m^2\n', params.gamma_lagb);
    fprintf('HAGB threshold angle: %.1f degrees\n', params.theta_HAGB);
    fprintf('===================================================\n\n');
end

function CompareAllMaterials()
    materials = {'LowCarbonSteel', 'AISI304', 'Aluminum6061', 'Copper', 'Titanium', 'Nickel718', 'Magnesium'};
    
    figure('Position', [100, 100, 1200, 800]);
    
    mu_vals = zeros(1, length(materials));
    Qact_vals = zeros(1, length(materials));
    Mgb_vals = zeros(1, length(materials));
    rho_c_vals = zeros(1, length(materials));
    
    for i = 1:length(materials)
        p = MaterialDatabase(materials{i});
        mu_vals(i) = p.mu / 1e9;
        Qact_vals(i) = p.Qact / 1000;
        Mgb_vals(i) = p.Mgb0;
        rho_c_vals(i) = p.rho_critical;
    end
    
    subplot(2,2,1);
    bar(mu_vals);
    set(gca, 'XTickLabel', materials, 'XTickLabelRotation', 45);
    ylabel('Shear Modulus (GPa)');
    title('Shear Modulus Comparison');
    grid on;
    
    subplot(2,2,2);
    bar(Qact_vals);
    set(gca, 'XTickLabel', materials, 'XTickLabelRotation', 45);
    ylabel('Activation Energy (kJ/mol)');
    title('Activation Energy Comparison');
    grid on;
    
    subplot(2,2,3);
    bar(log10(Mgb_vals));
    set(gca, 'XTickLabel', materials, 'XTickLabelRotation', 45);
    ylabel('log_{10}(M_{gb}) [m^4/(J*s)]');
    title('GB Mobility Prefactor Comparison');
    grid on;
    
    subplot(2,2,4);
    bar(rho_c_vals / 1e14);
    set(gca, 'XTickLabel', materials, 'XTickLabelRotation', 45);
    ylabel('Critical Dislocation Density (10^{14} m^{-2})');
    title('Critical Dislocation Density Comparison');
    grid on;
    
    sgtitle('Material Parameters Comparison', 'FontSize', 14);
end