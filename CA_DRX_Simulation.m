function CA_DRX_Simulation()
    clc; clear; close all;
    
    params = SetMaterialParameters();
    params = SetProcessParameters(params);
    params = SetCAParameters(params);
    
    [grainID, orientation, rho, isRecrystallized] = InitializeMicrostructure(params);
    
    results = struct();
    results.strain = [];
    results.stress = [];
    results.XdrxFraction = [];
    results.avgGrainSize = [];
    results.LAGBFraction = [];
    results.HAGBFraction = [];
    results.grainSizeDistribution = {};
    
    figure('Position', [100, 100, 1400, 900]);
    
    totalSteps = round(params.totalStrain / params.strainIncrement);
    currentStrain = 0;
    
    for step = 1:totalSteps
        currentStrain = currentStrain + params.strainIncrement;
        
        rho = UpdateDislocationDensity(rho, isRecrystallized, params, currentStrain);
        
        [grainID, orientation, rho, isRecrystallized] = Nucleation(grainID, orientation, rho, isRecrystallized, params);
        
        [grainID, orientation, rho, isRecrystallized] = GrainGrowth(grainID, orientation, rho, isRecrystallized, params);
        
        [grainID, orientation, rho, isRecrystallized] = TopologicalDeformation(grainID, orientation, rho, isRecrystallized, params, currentStrain);
        
        if mod(step, params.outputInterval) == 0 || step == totalSteps
            stats = CalculateStatistics(grainID, orientation, isRecrystallized, rho, params);
            stress = CalculateFlowStress(rho, params);
            
            results.strain(end+1) = currentStrain;
            results.stress(end+1) = stress;
            results.XdrxFraction(end+1) = stats.XdrxFraction;
            results.avgGrainSize(end+1) = stats.avgGrainSize;
            results.LAGBFraction(end+1) = stats.LAGBFraction;
            results.HAGBFraction(end+1) = stats.HAGBFraction;
            results.grainSizeDistribution{end+1} = stats.grainSizes;
            
            VisualizeMicrostructure(grainID, orientation, isRecrystallized, rho, results, params, currentStrain);
            drawnow;
        end
    end
    
    ExportResults(results, params);
    
    fprintf('\n===== Simulation Complete =====\n');
    fprintf('Final Strain: %.3f\n', currentStrain);
    fprintf('Final Stress: %.2f MPa\n', results.stress(end));
    fprintf('DRX Fraction: %.2f%%\n', results.XdrxFraction(end)*100);
    fprintf('Average Grain Size: %.2f um\n', results.avgGrainSize(end));
    fprintf('LAGB Fraction: %.2f%%\n', results.LAGBFraction(end)*100);
    fprintf('HAGB Fraction: %.2f%%\n', results.HAGBFraction(end)*100);
end

function params = SetMaterialParameters()
    params.materialName = 'Low Carbon Steel';
    
    params.b = 2.58e-10;
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
end

function params = SetProcessParameters(params)
    params.temperature = 1273;
    params.strainRate = 1.0;
    params.totalStrain = 0.8;
    
    params.strainIncrement = 0.001;
    params.outputInterval = 20;
end

function params = SetCAParameters(params)
    params.gridSize = 300;
    params.cellSize = 1e-6;
    
    params.initialGrainNum = 50;
    params.initialGrainSize = 50e-6;
    
    params.neighborhoodType = 'Moore';
    
    params.maxGrainID = 50000;
    
    params.deformationType = 'compression';
    params.deformationInterval = 100;
    
    if ~isfield(params, 'maxGrowthFraction')
        params.maxGrowthFraction = 0.02;
    end
    if ~isfield(params, 'rho_saturation')
        params.rho_saturation = 1e15;
    end
    if ~isfield(params, 'rho_critical')
        params.rho_critical = 1e13;
    end
end

function [grainID, orientation, rho, isRecrystallized] = InitializeMicrostructure(params)
    N = params.gridSize;
    grainID = zeros(N, N);
    orientation = zeros(N, N, 3);
    rho = ones(N, N) * params.rho0;
    isRecrystallized = false(N, N);
    
    numGrains = params.initialGrainNum;
    seeds_x = randi(N, numGrains, 1);
    seeds_y = randi(N, numGrains, 1);
    
    for i = 1:numGrains
        orientation(seeds_y(i), seeds_x(i), :) = rand(1, 3) * 360;
    end
    
    [X, Y] = meshgrid(1:N, 1:N);
    for i = 1:N
        for j = 1:N
            distances = sqrt((seeds_x - j).^2 + (seeds_y - i).^2);
            [~, closestGrain] = min(distances);
            grainID(i, j) = closestGrain;
            orientation(i, j, :) = orientation(seeds_y(closestGrain), seeds_x(closestGrain), :);
        end
    end
end

function rho = UpdateDislocationDensity(rho, isRecrystallized, params, strain)
    T = params.temperature;
    strainRate = params.strainRate;
    dEps = params.strainIncrement;
    
    Z = strainRate * exp(params.Qact / (params.R * T));
    k2 = params.k2_0 * (Z)^(-params.m_k2);
    
    drho = (params.k1 * sqrt(rho) - k2 * rho) * dEps;
    
    drho_unrex = drho;
    drho_rex = drho * 0.3;
    
    rho(~isRecrystallized) = rho(~isRecrystallized) + drho_unrex(~isRecrystallized);
    rho(isRecrystallized) = rho(isRecrystallized) + drho_rex(isRecrystallized);
    
    rho(rho < params.rho0) = params.rho0;
    rho(rho > params.rho_saturation) = params.rho_saturation;
end

function [grainID, orientation, rho, isRecrystallized] = Nucleation(grainID, orientation, rho, isRecrystallized, params)
    N = params.gridSize;
    T = params.temperature;
    
    [boundaryMask, ~] = FindGrainBoundaries(grainID);
    
    highRhoMask = rho > params.rho_critical;
    
    nucleationSites = boundaryMask & highRhoMask & ~isRecrystallized;
    
    baseProb = params.C_nuc * exp(-params.Qnuc / (params.R * T));
    
    [nucY, nucX] = find(nucleationSites);
    
    if isempty(nucY)
        return;
    end
    
    maxNewGrainID = max(grainID(:));
    
    idx = randperm(length(nucX));
    
    for k = 1:length(idx)
        i = nucY(idx(k));
        j = nucX(idx(k));
        
        if isRecrystallized(i, j)
            continue;
        end
        
        rhoRatio = rho(i,j) / params.rho_critical;
        nucProb = baseProb * params.nuc_probability_base * rhoRatio;
        nucProb = min(nucProb, 0.5);
        
        if rand() < nucProb
            maxNewGrainID = maxNewGrainID + 1;
            if maxNewGrainID > params.maxGrainID
                maxNewGrainID = params.maxGrainID;
            end
            
            grainID(i, j) = maxNewGrainID;
            orientation(i, j, :) = rand(1, 3) * 360;
            rho(i, j) = params.rho0;
            isRecrystallized(i, j) = true;
        end
    end
end

function [grainID, orientation, rho, isRecrystallized] = GrainGrowth(grainID, orientation, rho, isRecrystallized, params)
    N = params.gridSize;
    T = params.temperature;
    
    Mgb = params.Mgb0 * exp(-params.Qgb / (params.R * T));
    
    [boundaryMask, ~] = FindGrainBoundaries(grainID);
    
    rexGrainBoundary = boundaryMask & isRecrystallized;
    
    [boundaryY, boundaryX] = find(rexGrainBoundary);
    
    if isempty(boundaryY)
        return;
    end
    
    idx = randperm(length(boundaryX));
    
    di = [-1, -1, -1, 0, 0, 1, 1, 1];
    dj = [-1, 0, 1, -1, 1, -1, 0, 1];
    
    maxTransitions = round(N * N * params.maxGrowthFraction);
    transitionCount = 0;
    
    for k = 1:length(idx)
        if transitionCount >= maxTransitions
            break;
        end
        
        ii = boundaryY(idx(k));
        jj = boundaryX(idx(k));
        
        currentGrain = grainID(ii, jj);
        
        for n = 1:length(di)
            ni = ii + di(n);
            nj = jj + dj(n);
            
            if ni < 1 || ni > N || nj < 1 || nj > N
                continue;
            end
            
            if grainID(ni, nj) == currentGrain
                continue;
            end
            
            if isRecrystallized(ni, nj)
                continue;
            end
            
            neighborRho = rho(ni, nj);
            dRho = neighborRho - params.rho0;
            
            if dRho > params.rho_critical * 0.1
                dG = 0.5 * params.mu * params.b^2 * dRho;
                
                misori = CalculateMisorientation(squeeze(orientation(ii,jj,:)), squeeze(orientation(ni,nj,:)));
                if misori >= params.theta_HAGB
                    gamma = params.gamma_gb;
                else
                    gamma = params.gamma_lagb * (misori / params.theta_HAGB + 0.1);
                end
                
                P = dG - gamma / params.cellSize;
                
                if P > 0
                    velocity = Mgb * P;
                    dt = params.strainIncrement / params.strainRate;
                    transitionProb = velocity * dt / params.cellSize;
                    transitionProb = min(transitionProb, 0.8);
                    
                    if rand() < transitionProb
                        grainID(ni, nj) = currentGrain;
                        orientation(ni, nj, :) = orientation(ii, jj, :);
                        rho(ni, nj) = params.rho0;
                        isRecrystallized(ni, nj) = true;
                        transitionCount = transitionCount + 1;
                        break;
                    end
                end
            end
        end
    end
end

function [grainID, orientation, rho, isRecrystallized] = TopologicalDeformation(grainID, orientation, rho, isRecrystallized, params, currentStrain)
    persistent stepCounter
    if isempty(stepCounter)
        stepCounter = 0;
    end
    stepCounter = stepCounter + 1;
    
    if mod(stepCounter, params.deformationInterval) ~= 0
        return;
    end
    
    N = params.gridSize;
    
    switch params.deformationType
        case 'compression'
            strainFactor = params.strainIncrement * params.deformationInterval;
            
            newN_y = round(N * (1 - strainFactor));
            newN_x = round(N * (1 + strainFactor));
            
            if newN_y < N * 0.5 || newN_x > N * 2
                return;
            end
            
            [X_old, Y_old] = meshgrid(1:N, 1:N);
            [X_new, Y_new] = meshgrid(linspace(1, N, N), linspace(1, N, N));
            
            X_def = X_new * (N / newN_x) + (N - N^2/newN_x) / 2;
            Y_def = Y_new * (N / newN_y) + (N - N^2/newN_y) / 2;
            
            X_def = max(1, min(N, X_def));
            Y_def = max(1, min(N, Y_def));
            
            grainID_new = interp2(X_old, Y_old, double(grainID), X_def, Y_def, 'nearest');
            rho_new = interp2(X_old, Y_old, rho, X_def, Y_def, 'nearest');
            isRecrystallized_new = interp2(X_old, Y_old, double(isRecrystallized), X_def, Y_def, 'nearest') > 0.5;
            
            orientation_new = zeros(N, N, 3);
            for c = 1:3
                orientation_new(:,:,c) = interp2(X_old, Y_old, orientation(:,:,c), X_def, Y_def, 'nearest');
            end
            
            grainID = round(grainID_new);
            orientation = orientation_new;
            rho = rho_new;
            isRecrystallized = isRecrystallized_new;
            
        case 'tension'
            strainFactor = params.strainIncrement * params.deformationInterval;
            
            [X_old, Y_old] = meshgrid(1:N, 1:N);
            [X_new, Y_new] = meshgrid(linspace(1, N, N), linspace(1, N, N));
            
            X_def = X_new * (1 + strainFactor);
            Y_def = Y_new * (1 - strainFactor * 0.5);
            
            X_def = max(1, min(N, X_def));
            Y_def = max(1, min(N, Y_def));
            
            grainID = round(interp2(X_old, Y_old, double(grainID), X_def, Y_def, 'nearest'));
            rho = interp2(X_old, Y_old, rho, X_def, Y_def, 'nearest');
            isRecrystallized = interp2(X_old, Y_old, double(isRecrystallized), X_def, Y_def, 'nearest') > 0.5;
            
            for c = 1:3
                orientation(:,:,c) = interp2(X_old, Y_old, orientation(:,:,c), X_def, Y_def, 'nearest');
            end
            
        case 'shear'
            strainFactor = params.strainIncrement * params.deformationInterval;
            
            [X_old, Y_old] = meshgrid(1:N, 1:N);
            
            X_def = X_old + strainFactor * (Y_old - N/2);
            Y_def = Y_old;
            
            X_def = mod(X_def - 1, N) + 1;
            
            grainID = round(interp2(X_old, Y_old, double(grainID), X_def, Y_def, 'nearest'));
            rho = interp2(X_old, Y_old, rho, X_def, Y_def, 'nearest');
            isRecrystallized = interp2(X_old, Y_old, double(isRecrystallized), X_def, Y_def, 'nearest') > 0.5;
            
            for c = 1:3
                orientation(:,:,c) = interp2(X_old, Y_old, orientation(:,:,c), X_def, Y_def, 'nearest');
            end
    end
    
    grainID(isnan(grainID)) = 1;
    rho(isnan(rho)) = params.rho0;
end

function [boundaryMask, boundaryType] = FindGrainBoundaries(grainID)
    N = size(grainID, 1);
    boundaryMask = false(N, N);
    boundaryType = zeros(N, N);
    
    for i = 2:N-1
        for j = 2:N-1
            neighbors = [grainID(i-1,j), grainID(i+1,j), grainID(i,j-1), grainID(i,j+1)];
            if any(neighbors ~= grainID(i,j))
                boundaryMask(i,j) = true;
            end
        end
    end
end

function neighbors = GetNeighbors(i, j, N, type)
    switch type
        case 'VonNeumann'
            di = [-1, 1, 0, 0];
            dj = [0, 0, -1, 1];
        case 'Moore'
            di = [-1, -1, -1, 0, 0, 1, 1, 1];
            dj = [-1, 0, 1, -1, 1, -1, 0, 1];
    end
    
    neighbors = [];
    for k = 1:length(di)
        ni = i + di(k);
        nj = j + dj(k);
        if ni >= 1 && ni <= N && nj >= 1 && nj <= N
            neighbors = [neighbors; ni, nj];
        end
    end
end

function misorientation = CalculateMisorientation(ori1, ori2)
    dOri = abs(ori1 - ori2);
    dOri = min(dOri, 360 - dOri);
    misorientation = sqrt(sum(dOri.^2)) / sqrt(3);
end

function stress = CalculateFlowStress(rho, params)
    avgRho = mean(rho(:));
    avgRho = min(avgRho, params.rho_saturation);
    stress = params.alpha_taylor * params.mu * params.b * sqrt(avgRho);
    stress = stress / 1e6;
    
    T = params.temperature;
    Tm = 1800;
    tempFactor = 1 - 0.5 * (T / Tm);
    tempFactor = max(0.3, tempFactor);
    stress = stress * tempFactor;
end

function stats = CalculateStatistics(grainID, orientation, isRecrystallized, rho, params)
    N = params.gridSize;
    
    stats.XdrxFraction = sum(isRecrystallized(:)) / numel(isRecrystallized);
    
    uniqueGrains = unique(grainID(:));
    grainSizes = zeros(length(uniqueGrains), 1);
    grainRecrystallized = false(length(uniqueGrains), 1);
    
    for k = 1:length(uniqueGrains)
        mask = grainID == uniqueGrains(k);
        grainSizes(k) = sqrt(sum(mask(:)) * params.cellSize^2 * 4 / pi) * 1e6;
        grainRecrystallized(k) = any(isRecrystallized(mask));
    end
    
    stats.grainSizes = grainSizes;
    stats.avgGrainSize = mean(grainSizes);
    stats.stdGrainSize = std(grainSizes);
    stats.drxGrainSize = mean(grainSizes(grainRecrystallized));
    stats.unrexGrainSize = mean(grainSizes(~grainRecrystallized));
    
    if isnan(stats.drxGrainSize), stats.drxGrainSize = 0; end
    if isnan(stats.unrexGrainSize), stats.unrexGrainSize = stats.avgGrainSize; end
    
    [boundaryMask, ~] = FindGrainBoundaries(grainID);
    [by, bx] = find(boundaryMask);
    
    LAGBCount = 0;
    HAGBCount = 0;
    totalBoundary = 0;
    
    for k = 1:length(bx)
        i = by(k);
        j = bx(k);
        
        neighbors = GetNeighbors(i, j, N, 'VonNeumann');
        
        for n = 1:size(neighbors, 1)
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            
            if grainID(ni, nj) ~= grainID(i, j)
                misori = CalculateMisorientation(squeeze(orientation(i,j,:)), squeeze(orientation(ni,nj,:)));
                totalBoundary = totalBoundary + 1;
                
                if misori >= params.theta_HAGB
                    HAGBCount = HAGBCount + 1;
                else
                    LAGBCount = LAGBCount + 1;
                end
            end
        end
    end
    
    if totalBoundary > 0
        stats.LAGBFraction = LAGBCount / totalBoundary;
        stats.HAGBFraction = HAGBCount / totalBoundary;
    else
        stats.LAGBFraction = 0;
        stats.HAGBFraction = 1;
    end
    
    stats.totalGrains = length(uniqueGrains);
    stats.avgDislocationDensity = mean(rho(:));
end

function VisualizeMicrostructure(grainID, orientation, isRecrystallized, rho, results, params, currentStrain)
    subplot(2,3,1);
    imagesc(grainID);
    colormap(gca, jet(256));
    axis equal tight;
    title(sprintf('Grain Structure (\\epsilon = %.3f)', currentStrain));
    colorbar;
    
    subplot(2,3,2);
    IPFColor = zeros(size(grainID, 1), size(grainID, 2), 3);
    IPFColor(:,:,1) = orientation(:,:,1) / 360;
    IPFColor(:,:,2) = orientation(:,:,2) / 360;
    IPFColor(:,:,3) = orientation(:,:,3) / 360;
    imagesc(IPFColor);
    axis equal tight;
    title('Orientation (IPF-like)');
    
    subplot(2,3,3);
    rhoDisplay = log10(rho);
    imagesc(rhoDisplay);
    colormap(gca, hot);
    axis equal tight;
    title('Dislocation Density (log_{10})');
    colorbar;
    
    subplot(2,3,4);
    if ~isempty(results.strain)
        yyaxis left;
        plot(results.strain, results.stress, 'b-', 'LineWidth', 2);
        ylabel('Flow Stress (MPa)');
        yyaxis right;
        plot(results.strain, results.XdrxFraction * 100, 'r-', 'LineWidth', 2);
        ylabel('DRX Fraction (%)');
        xlabel('Strain');
        title('Stress-Strain & DRX Kinetics');
        legend('Stress', 'X_{DRX}', 'Location', 'best');
        grid on;
    end
    
    subplot(2,3,5);
    if ~isempty(results.grainSizeDistribution) && ~isempty(results.grainSizeDistribution{end})
        histogram(results.grainSizeDistribution{end}, 20, 'FaceColor', [0.3 0.6 0.9]);
        xlabel('Grain Size (\mum)');
        ylabel('Count');
        title(sprintf('Grain Size Distribution (Avg: %.1f \\mum)', results.avgGrainSize(end)));
        grid on;
    end
    
    subplot(2,3,6);
    if ~isempty(results.LAGBFraction)
        bar([results.LAGBFraction(end), results.HAGBFraction(end)] * 100);
        set(gca, 'XTickLabel', {'LAGB (<15°)', 'HAGB (≥15°)'});
        ylabel('Fraction (%)');
        title('Grain Boundary Character');
        ylim([0 100]);
        grid on;
    end
    
    sgtitle(sprintf('%s: T = %d K, \\epsilon_{rate} = %.2f s^{-1}', ...
        params.materialName, params.temperature, params.strainRate), 'FontSize', 12);
end

function ExportResults(results, params)
    T = table(results.strain', results.stress', results.XdrxFraction'*100, ...
        results.avgGrainSize', results.LAGBFraction'*100, results.HAGBFraction'*100, ...
        'VariableNames', {'Strain', 'Stress_MPa', 'DRX_Fraction_pct', ...
        'AvgGrainSize_um', 'LAGB_pct', 'HAGB_pct'});
    
    filename = sprintf('DRX_Results_T%d_SR%.2f.csv', params.temperature, params.strainRate);
    writetable(T, filename);
    fprintf('Results exported to %s\n', filename);
    
    save('DRX_Simulation_Results.mat', 'results', 'params');
    fprintf('Full results saved to DRX_Simulation_Results.mat\n');
end