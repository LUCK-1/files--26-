function ParameterSweep_DRX()
    clc; clear; close all;
    
    temperatures = [1173, 1273, 1373];
    strainRates = [0.1, 1.0, 10.0];
    totalStrains = [0.6, 0.8, 1.0];
    
    allResults = {};
    resultIndex = 1;
    
    for T_idx = 1:length(temperatures)
        for SR_idx = 1:length(strainRates)
            for S_idx = 1:length(totalStrains)
                T = temperatures(T_idx);
                SR = strainRates(SR_idx);
                S = totalStrains(S_idx);
                
                fprintf('\n========================================\n');
                fprintf('Running: T=%dK, SR=%.1f/s, Strain=%.1f\n', T, SR, S);
                fprintf('========================================\n');
                
                params = SetMaterialParameters();
                params.temperature = T;
                params.strainRate = SR;
                params.totalStrain = S;
                params = SetCAParameters(params);
                
                result = RunSingleSimulation(params);
                
                allResults{resultIndex} = result;
                allResults{resultIndex}.temperature = T;
                allResults{resultIndex}.strainRate = SR;
                allResults{resultIndex}.totalStrain = S;
                resultIndex = resultIndex + 1;
            end
        end
    end
    
    PlotParameterStudy(allResults, temperatures, strainRates);
    
    save('ParameterSweep_Results.mat', 'allResults', 'temperatures', 'strainRates', 'totalStrains');
    fprintf('\nAll results saved to ParameterSweep_Results.mat\n');
end

function params = SetMaterialParameters()
    params.materialName = 'Low Carbon Steel';
    params.b = 2.58e-10;
    params.mu = 80e9;
    params.Qact = 280e3;
    params.R = 8.314;
    params.k1 = 2.5e9;
    params.k2_0 = 5;
    params.m_k2 = 0.12;
    params.Qk2 = 120e3;
    params.rho0 = 1e12;
    params.rho_critical = 8e13;
    params.Qnuc = 150e3;
    params.C_nuc = 5e-2;
    params.nuc_probability_base = 0.15;
    params.Qgb = 130e3;
    params.Mgb0 = 5e-5;
    params.gamma_gb = 0.5;
    params.gamma_lagb = 0.15;
    params.theta_HAGB = 15;
    params.alpha_taylor = 0.5;
    params.strainIncrement = 0.001;
    params.outputInterval = 50;
end

function params = SetCAParameters(params)
    params.gridSize = 200;
    params.cellSize = 1e-6;
    params.initialGrainNum = 30;
    params.initialGrainSize = 50e-6;
    params.neighborhoodType = 'Moore';
    params.maxGrainID = 10000;
    params.deformationType = 'compression';
    params.deformationInterval = 50;
end

function result = RunSingleSimulation(params)
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
    
    result = struct();
    result.strain = [];
    result.stress = [];
    result.XdrxFraction = [];
    result.avgGrainSize = [];
    result.LAGBFraction = [];
    result.HAGBFraction = [];
    
    totalSteps = round(params.totalStrain / params.strainIncrement);
    currentStrain = 0;
    
    for step = 1:totalSteps
        currentStrain = currentStrain + params.strainIncrement;
        
        T = params.temperature;
        strainRate = params.strainRate;
        dEps = params.strainIncrement;
        k2 = params.k2_0 * (strainRate * exp(params.Qk2 / (params.R * T)))^params.m_k2;
        drho = (params.k1 * sqrt(rho) - k2 * rho) * dEps;
        rho = rho + drho;
        rho(rho < params.rho0) = params.rho0;
        rho(isRecrystallized) = rho(isRecrystallized) * 0.98;
        rho(isRecrystallized & rho < params.rho0) = params.rho0;
        
        boundaryMask = false(N, N);
        for i = 2:N-1
            for j = 2:N-1
                neighbors = [grainID(i-1,j), grainID(i+1,j), grainID(i,j-1), grainID(i,j+1)];
                if any(neighbors ~= grainID(i,j))
                    boundaryMask(i,j) = true;
                end
            end
        end
        
        highRhoMask = rho > params.rho_critical;
        nucleationSites = boundaryMask & highRhoMask & ~isRecrystallized;
        nucleationProb = params.C_nuc * exp(-params.Qnuc / (params.R * T));
        nucleationProb = nucleationProb * (rho / params.rho_critical);
        
        [nucY, nucX] = find(nucleationSites);
        maxNewGrainID = max(grainID(:));
        
        for k = 1:length(nucX)
            ii = nucY(k);
            jj = nucX(k);
            if rand() < nucleationProb(ii, jj) * params.nuc_probability_base
                maxNewGrainID = maxNewGrainID + 1;
                grainID(ii, jj) = maxNewGrainID;
                orientation(ii, jj, :) = rand(1, 3) * 360;
                rho(ii, jj) = params.rho0;
                isRecrystallized(ii, jj) = true;
            end
        end
        
        Mgb = params.Mgb0 * exp(-params.Qgb / (params.R * T));
        [boundaryY, boundaryX] = find(boundaryMask);
        idx = randperm(length(boundaryX));
        
        for k = 1:length(idx)
            ii = boundaryY(idx(k));
            jj = boundaryX(idx(k));
            
            di = [-1, -1, -1, 0, 0, 1, 1, 1];
            dj = [-1, 0, 1, -1, 1, -1, 0, 1];
            
            for nn = 1:length(di)
                ni = ii + di(nn);
                nj = jj + dj(nn);
                if ni >= 1 && ni <= N && nj >= 1 && nj <= N
                    if grainID(ni, nj) ~= grainID(ii, jj)
                        dRho = rho(ii, jj) - rho(ni, nj);
                        if dRho > 0
                            dG = 0.5 * params.mu * params.b^2 * dRho;
                            ori1 = squeeze(orientation(ii,jj,:));
                            ori2 = squeeze(orientation(ni,nj,:));
                            dOri = abs(ori1 - ori2);
                            dOri = min(dOri, 360 - dOri);
                            misori = sqrt(sum(dOri.^2)) / sqrt(3);
                            
                            if misori >= params.theta_HAGB
                                gamma = params.gamma_gb;
                            else
                                gamma = params.gamma_lagb * (misori / params.theta_HAGB);
                            end
                            
                            P = dG - gamma / params.cellSize;
                            if P > 0
                                velocity = Mgb * P;
                                transProb = velocity * params.strainIncrement / params.strainRate / params.cellSize;
                                transProb = min(transProb, 0.5);
                                if rand() < transProb
                                    grainID(ii, jj) = grainID(ni, nj);
                                    orientation(ii, jj, :) = orientation(ni, nj, :);
                                    rho(ii, jj) = rho(ni, nj);
                                    isRecrystallized(ii, jj) = isRecrystallized(ni, nj);
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if mod(step, params.outputInterval) == 0 || step == totalSteps
            XdrxFrac = sum(isRecrystallized(:)) / numel(isRecrystallized);
            avgRho = mean(rho(:));
            stress = params.alpha_taylor * params.mu * params.b * sqrt(avgRho) / 1e6;
            
            uniqueGrains = unique(grainID(:));
            grainSizes = zeros(length(uniqueGrains), 1);
            for kk = 1:length(uniqueGrains)
                mask = grainID == uniqueGrains(kk);
                grainSizes(kk) = sqrt(sum(mask(:)) * params.cellSize^2 * 4 / pi) * 1e6;
            end
            avgGS = mean(grainSizes);
            
            LAGBCount = 0;
            HAGBCount = 0;
            totalBnd = 0;
            [by, bx] = find(boundaryMask);
            for kk = 1:min(length(bx), 1000)
                ii = by(kk);
                jj = bx(kk);
                for dd = 1:4
                    if dd == 1, ni = ii-1; nj = jj; end
                    if dd == 2, ni = ii+1; nj = jj; end
                    if dd == 3, ni = ii; nj = jj-1; end
                    if dd == 4, ni = ii; nj = jj+1; end
                    if ni >= 1 && ni <= N && nj >= 1 && nj <= N
                        if grainID(ni, nj) ~= grainID(ii, jj)
                            ori1 = squeeze(orientation(ii,jj,:));
                            ori2 = squeeze(orientation(ni,nj,:));
                            dOri = abs(ori1 - ori2);
                            dOri = min(dOri, 360 - dOri);
                            misori = sqrt(sum(dOri.^2)) / sqrt(3);
                            totalBnd = totalBnd + 1;
                            if misori >= params.theta_HAGB
                                HAGBCount = HAGBCount + 1;
                            else
                                LAGBCount = LAGBCount + 1;
                            end
                        end
                    end
                end
            end
            
            if totalBnd > 0
                LAGBFrac = LAGBCount / totalBnd;
                HAGBFrac = HAGBCount / totalBnd;
            else
                LAGBFrac = 0;
                HAGBFrac = 1;
            end
            
            result.strain(end+1) = currentStrain;
            result.stress(end+1) = stress;
            result.XdrxFraction(end+1) = XdrxFrac;
            result.avgGrainSize(end+1) = avgGS;
            result.LAGBFraction(end+1) = LAGBFrac;
            result.HAGBFraction(end+1) = HAGBFrac;
        end
    end
    
    result.finalGrainID = grainID;
    result.finalOrientation = orientation;
    result.finalRho = rho;
    result.finalIsRex = isRecrystallized;
end

function PlotParameterStudy(allResults, temperatures, strainRates)
    figure('Position', [50, 50, 1600, 1000]);
    
    colors = lines(length(temperatures));
    lineStyles = {'-', '--', ':'};
    
    subplot(2,3,1);
    hold on;
    for i = 1:length(allResults)
        T_idx = find(temperatures == allResults{i}.temperature);
        SR_idx = find(strainRates == allResults{i}.strainRate);
        plot(allResults{i}.strain, allResults{i}.stress, ...
            'Color', colors(T_idx,:), 'LineStyle', lineStyles{SR_idx}, 'LineWidth', 1.5);
    end
    xlabel('Strain');
    ylabel('Flow Stress (MPa)');
    title('Flow Stress vs Strain');
    grid on;
    
    subplot(2,3,2);
    hold on;
    for i = 1:length(allResults)
        T_idx = find(temperatures == allResults{i}.temperature);
        SR_idx = find(strainRates == allResults{i}.strainRate);
        plot(allResults{i}.strain, allResults{i}.XdrxFraction*100, ...
            'Color', colors(T_idx,:), 'LineStyle', lineStyles{SR_idx}, 'LineWidth', 1.5);
    end
    xlabel('Strain');
    ylabel('DRX Fraction (%)');
    title('DRX Kinetics');
    grid on;
    
    subplot(2,3,3);
    hold on;
    for i = 1:length(allResults)
        T_idx = find(temperatures == allResults{i}.temperature);
        SR_idx = find(strainRates == allResults{i}.strainRate);
        plot(allResults{i}.strain, allResults{i}.avgGrainSize, ...
            'Color', colors(T_idx,:), 'LineStyle', lineStyles{SR_idx}, 'LineWidth', 1.5);
    end
    xlabel('Strain');
    ylabel('Average Grain Size (\mum)');
    title('Grain Size Evolution');
    grid on;
    
    subplot(2,3,4);
    finalXdrx = zeros(length(temperatures), length(strainRates));
    for i = 1:length(allResults)
        T_idx = find(temperatures == allResults{i}.temperature);
        SR_idx = find(strainRates == allResults{i}.strainRate);
        finalXdrx(T_idx, SR_idx) = allResults{i}.XdrxFraction(end) * 100;
    end
    bar3(finalXdrx);
    set(gca, 'XTickLabel', strainRates);
    set(gca, 'YTickLabel', temperatures);
    xlabel('Strain Rate (s^{-1})');
    ylabel('Temperature (K)');
    zlabel('Final DRX Fraction (%)');
    title('Parameter Effect on DRX');
    
    subplot(2,3,5);
    finalGS = zeros(length(temperatures), length(strainRates));
    for i = 1:length(allResults)
        T_idx = find(temperatures == allResults{i}.temperature);
        SR_idx = find(strainRates == allResults{i}.strainRate);
        finalGS(T_idx, SR_idx) = allResults{i}.avgGrainSize(end);
    end
    bar3(finalGS);
    set(gca, 'XTickLabel', strainRates);
    set(gca, 'YTickLabel', temperatures);
    xlabel('Strain Rate (s^{-1})');
    ylabel('Temperature (K)');
    zlabel('Final Grain Size (\mum)');
    title('Parameter Effect on Grain Size');
    
    subplot(2,3,6);
    hold on;
    for i = 1:length(allResults)
        T_idx = find(temperatures == allResults{i}.temperature);
        SR_idx = find(strainRates == allResults{i}.strainRate);
        plot(allResults{i}.strain, allResults{i}.HAGBFraction*100, ...
            'Color', colors(T_idx,:), 'LineStyle', lineStyles{SR_idx}, 'LineWidth', 1.5);
    end
    xlabel('Strain');
    ylabel('HAGB Fraction (%)');
    title('HAGB Evolution');
    grid on;
    
    legendStr = {};
    for T = temperatures
        for SR = strainRates
            legendStr{end+1} = sprintf('T=%dK, SR=%.1f', T, SR);
        end
    end
    legend(legendStr, 'Location', 'bestoutside', 'FontSize', 7);
    
    sgtitle('Parameter Study: Temperature and Strain Rate Effects on DRX', 'FontSize', 14);
    
    saveas(gcf, 'ParameterStudy_Results.png');
    fprintf('Figure saved to ParameterStudy_Results.png\n');
end
