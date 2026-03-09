%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %%%     %%%    %%%%%%%   %%%     %
%                                      %%%     %%%   %%%        %%%     %
%                                      %%%     %%%   %%%        %%%     %
%                                      %%%     %%%   %%%        %%%     %
%                                       %%%%%%%%%     %%%%%%%   %%%%%%%%%
%                                                                       %
%                                                                       %
%                      HypoSubstrateTesting                             %
% Testing hypothesized neural substrates with lesions                   %
%                                                                       %
% This script contains the function for plotting the distribution of    %
% p-values from simulation runs using split violin plots.               %
%                                                                       %
%                                                                       %
% last update: 2026/03/07.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_simulation_results_violin(varargin)
    % plot_simulation_results_violin - Create figures from simulation results
    % Uses custom split violin plots for binned p-value comparison
    %
    % Usage:
    %   plot_simulation_results_violin(resultFile)
    %   plot_simulation_results_violin(resultFile, 'Name', Value, ...)
    %
    % Arguments:
    %   resultFile  - Path to simulation results .mat file (from simulations_clean.m)
    %
    % Options:
    %   'fontPath'    - Path to DINPro font directory or font name (default: 'DINPro-Regular')
    %   'outputFile'  - Output filepath for PDF (default: [resultFile(1:end-4) '_figures.pdf'])
    %   'dataDir'     - Data directory containing required template/morphology files (default: '')
    %   'binSize'     - Noise bin width for split violin plot (default: from .mat file, or 0.05)
    %
    % Example:
    %   plot_simulation_results_violin('base-MDN_morphology.mat')
    %   plot_simulation_results_violin('base-MDN_morphology.mat', 'outputFile', 'output/figures.pdf')
    %   plot_simulation_results_violin('base-MDN_morphology.mat', 'binSize', 0.1)

    p = inputParser;
    %addParameter(p, 'fontPath', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'outputFile', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'dataDir', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'binSize', NaN, @(x) isnumeric(x) && isscalar(x));
    
    if nargin == 1
        resultFile = varargin{1};
        %fontPath = '';
        outputFile = '';
        dataDir = '';
        binSizeParam = NaN;
    else
        resultFile = varargin{1};
        parse(p, varargin{2:end});
        %fontPath = char(p.Results.fontPath);
        outputFile = char(p.Results.outputFile);
        dataDir = char(p.Results.dataDir);
        binSizeParam = p.Results.binSize;
    end
    
    if isstring(resultFile)
        resultFile = char(resultFile);
    end
    
    % Default font path
    %if isempty(fontPath)
    %    fontPath = 'H:\GitHub\hypolesions\dinprottf';
    %end
    
    fontName = "DINPro-Regular";

    fprintf('Loading simulation results from: %s\n', resultFile);
    load(resultFile);

    % Resolve binSize: parameter override > .mat file > default 0.05
    if ~isnan(binSizeParam)
        binSize = binSizeParam;
    elseif ~exist('binSize', 'var')
        binSize = 0.05;
    end

    fprintf('Results loaded:\n');
    fprintf('  Base image: %s\n', baseImageFile);
    fprintf('  Target: %s\n', regressionTarget);
    fprintf('  Iterations: %d\n', nits);
    fprintf('  Cases per iteration: %d\n', ncas);
    % fprintf('  Nnn-zero ratio: %d\n', nonZeroRatio)
    fprintf('  Noise range: %.2f - %.2f\n', min(nlev), max(nlev));

    [~, baseName, ~] = fileparts(baseImageFile);
    baseNameDisplay = strrep(baseName, '_', ' ');
    baseNameDisplay = regexprep(baseNameDisplay, '^\d+', '');
    baseNameDisplay = strrep(baseNameDisplay, 'MNI152', '');
    baseNameDisplay = strtrim(baseNameDisplay);

    fig = figure('Name', ['Simulation Results - ' upper(regressionTarget) ' (Base: ' baseNameDisplay ')'], ...
        'Visible', 'on', 'Position', [100 100 1600 1200]);
    set(fig, 'PaperUnits', 'points', 'PaperSize', [1600 1200], 'PaperPosition', [0 0 1600 1200]);
    tiledlayout('flow', 'TileSpacing', 'tight');
    set(gca, 'Color', 'white');

    nlev_plot = [min(nlev) max(nlev)];

    % Use res and princ from result file if available (saved by simulations.m)
    % Fall back to loading from separate files for backward compatibility
    if ~exist('res', 'var') || isempty(res)
        fprintf('res not found in result file, loading from template file...\n');
        res = computeResiduals(baseImageFile, regressionTarget, targetFile, dataDir);
    else
        fprintf('Using res from result file (%dx%d)\n', size(res, 1), size(res, 2));
    end
    nCases = size(res, 1);
    
    if ~exist('princ', 'var') || isempty(princ)
        fprintf('princ not found in result file, loading from source file...\n');
        princ = computePrinc(regressionTarget, targetFile, baseName, dataDir);
    else
        fprintf('Using princ from result file (%dx%d)\n', size(princ, 1), size(princ, 2));
    end

    dataValid = false;
    if ~isempty(res) && ~isempty(princ)
        if size(res, 1) == length(princ(:))
            dataValid = true;
        end
    end

    if dataValid
        for j = 1:length(nlev_plot)
            nexttile
            ax1 = gca;
            set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');

            noiseLevel = nlev_plot(j);

            selvec = zeros(nCases, 1);
            selvec(1:ncas) = 1;
            selvec = selvec(randperm(length(selvec)));
            indvec = logical([selvec; selvec]);

            BaseVol = res(:,1);
            nBaseVol = res(:,2);
            totvol = res(:,3);

            pd = makedist('Normal', 'mu', 0, 'sigma', noiseLevel);
            dis = random(pd, [nCases 1]);

            switch lower(regressionTarget)
                case 'morphology'
                    cog = 1 - rescale(princ(:));
                    cog = cog(:) + dis;
                case 'targetfile'
                    princ_temp = rescale(princ(:));
                    cog = 1 - princ_temp;
                    cog = cog(:) + dis;
                case 'baseimage'
                    princ_temp = rescale(res(:,1));
                    cog = 1 - princ_temp;
                    cog = cog + dis;
                    cog = cog(:);
            end

            cogBase = cog;
            totvolBase = totvol;
            cognBase = cog;
            totvolnBase = totvol;

            ovrl = [BaseVol; nBaseVol];
            voll = [totvolBase; totvolnBase];
            cogg = [cogBase; cognBase];
            comp = zeros(size(cogg, 1), 1);
            comp(1:nCases) = 1;

            data = table(cogg, ovrl, voll, comp, 'VariableNames', {'Cognition', 'Overlap', 'Volume', 'Network'});
            data.comp = categorical(data.Network);
            data(~indvec, :) = [];

            fit = fitlm(data, 'interactions', 'ResponseVar', 'Cognition', ...
                'PredictorVars', {'Overlap', 'Network'}, ...
                'CategoricalVar', {'Network'}, 'RobustOpts', 'off');

            w = linspace(min(ovrl), max(voll));

            gscatter(ovrl(indvec), cogg(indvec), comp(indvec), [0 0 0; 1 0 0], '..')
            line(w, feval(fit, w, 1), 'Color', [1 0 0], 'LineWidth', 2)
            line(w, feval(fit, w, 0), 'Color', 'k', 'LineWidth', 2)

            xlim([0 max(voll)]);
            xlabel('Compartment damage volume')
            ylabel('Simulated cognitive performance')
            title(['Noise \sigma = ' num2str(nlev_plot(j))]);
            axis square
            set(gca, 'XScale', 'linear', 'YScale', 'linear')

            anv = anova(fit);
            tmp = anv.pValue;

            text(100, 0, ['p=' num2str(tmp(3))], 'FontSize', 14, 'FontName', fontName)
            legend({['non-' baseNameDisplay ' overlap'], [baseNameDisplay ' overlap'], [baseNameDisplay ' fit'], ['non-' baseNameDisplay ' fit']}, ...
                'FontSize', 14, 'Box', 'off');
            set(gca, 'FontName', fontName, 'FontSize', 14)
            axis square
            ylim([-2.5 2.5]);
            xlim([0 300]);
            set(gcf, 'Color', 'w')
            box on;
        end
    else
        fprintf('Warning: Data mismatch (res: %d rows, princ: %d elements). Skipping scatter plots.\n', ...
            size(res, 1), length(princ(:)));
        nexttile; hold on;
        text(0.5, 0.5, 'Scatter plots skipped: data mismatch', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
        axis on; box on;
        nexttile; hold on;
        text(0.5, 0.5, 'Scatter plots skipped: data mismatch', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
        axis on; box on;
    end

    set(0, 'DefaultFigureColor', 'white')

    switch lower(regressionTarget)
        case 'morphology'
            figureTitle = ['Spurious ' baseNameDisplay ' association vs noise'];
        case 'targetfile'
            figureTitle = ['Spurious ' baseNameDisplay ' association vs noise'];
        otherwise
            figureTitle = [baseNameDisplay ' association vs noise'];
    end

    nexttile
    ax1 = gca;
    set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    box on;
    set(gca, 'Color', 'white');
    s1 = scatter(nlev, result(1, :), 10, [0 0.447 0.741], 'x', 'LineWidth', 0.5);
    hold on
    s2 = scatter(nlev, result(2, :), 10, [0.466 0.674 0.188], 'filled', 'MarkerEdgeColor', 'none');
    hold on
    title(figureTitle);
    xlabel('Noise (\sigma)')
    ylabel('-log(p)')
    xticks(min(nlev):0.05:max(nlev))
    xtickangle(0)
    plot([min(nlev) max(nlev)], [-log(0.05) -log(0.05)], 'k:', 'LineWidth', 2)
    legend({'Slope comparison', ['Slope comparison with' newline 'lesion volume covariate'], 'Significance threshold'}, ...
        'FontSize', 14, 'Box', 'off');
    set(gca, 'FontName', fontName, 'FontSize', 14)
    ylim([0 max([25, max(result(:)), max(result(:))*1.1])])
    axis square

    nexttile
    ax2 = gca;
    set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    box on;
    set(gca, 'Color', 'white');
    hold on

    % --- Compute noise bins (same logic as simulations.m CSV export) ---
    nlow = min(nlev);
    nhigh = max(nlev);
    binEdges = round(nlow:binSize:nhigh, 10);
    if binEdges(end) < nhigh
        binEdges(end+1) = round(binEdges(end) + binSize, 10);
    end
    nBins = length(binEdges) - 1;
    binIdx = discretize(nlev, binEdges, 'IncludedEdge', 'right');

    significanceThreshold = -log(0.05);
    colorLeft  = [0 0.447 0.741];     % blue  - without lesion adjustment
    colorRight = [0.466 0.674 0.188]; % green - with lesion adjustment
    faceAlpha  = 0.3;
    maxHalfWidth = 0.4;  % max half-width in x-axis units

    % --- Draw split violins for each bin ---
    binCenters = zeros(1, nBins);
    binLabels  = cell(1, nBins);
    hLeft  = [];
    hRight = [];
    pctLeftAll  = zeros(1, nBins);
    pctRightAll = zeros(1, nBins);
    q90All      = zeros(1, nBins);

    for b = 1:nBins
        inBin = (binIdx == b);
        dataLeft  = result(1, inBin);   % without lesion adjustment
        dataRight = result(2, inBin);   % with lesion adjustment

        binCenters(b) = b;
        binLabels{b} = sprintf('%.2f-%.2f', binEdges(b), binEdges(b+1));

        hL = drawHalfViolin(b, dataLeft, 'left', colorLeft, colorLeft, faceAlpha, maxHalfWidth);
        hR = drawHalfViolin(b, dataRight, 'right', colorRight, colorRight, faceAlpha, maxHalfWidth);

        if isempty(hLeft)
            hLeft = hL;
            hRight = hR;
        end

        % Compute per-bin stats for text annotations
        nInBin = sum(inBin);
        pctLeftAll(b)  = 100 * sum(dataLeft > significanceThreshold) / nInBin;
        pctRightAll(b) = 100 * sum(dataRight > significanceThreshold) / nInBin;
        q90All(b) = quantile(dataLeft, 0.9);
    end

    % Place all text labels at a uniform height: max of the 0.9 quantiles
    textY = max(q90All);
    for b = 1:nBins
        text(b - 0.15, textY, sprintf('%.0f%%', pctLeftAll(b)), ...
            'FontSize', 14, 'FontName', fontName, 'Color', colorLeft, ...
            'HorizontalAlignment', 'right');
        text(b + 0.15, textY, sprintf('%.0f%%', pctRightAll(b)), ...
            'FontSize', 14, 'FontName', fontName, 'Color', colorRight, ...
            'HorizontalAlignment', 'left');
    end

    % Significance threshold line
    hThresh = plot([0.5 nBins+0.5], [significanceThreshold significanceThreshold], ...
        'k:', 'LineWidth', 2, 'DisplayName', 'Significance threshold');

    xticks(binCenters);
    xticklabels(binLabels);
    xtickangle(45)
    xlim([0.5 nBins + 0.5])
    xlabel('Noise bin (\sigma)')
    ylabel('-log(p)')
    title('Effect of lesion volume adjustment');

    % Legend using representative handles
    legend([hLeft, hRight, hThresh], ...
        {'Slope comparison', 'with lesion volume covariate', 'Significance threshold'}, ...
        'FontSize', 14, 'Box', 'off', 'Location', 'northeast');
    set(gca, 'FontName', fontName, 'FontSize', 14)
    ylim([0 max([25, max(result(:)), max(result(:))*1.1])])
    set(gcf, 'Color', 'w')

    hold off;

    if isempty(outputFile)
        [~, baseNameOut, ~] = fileparts(resultFile);
        outputFile = [baseNameOut '_figures.pdf'];
    end

    fprintf('\nSaving figure to: %s\n', outputFile);
    [outPath, ~, ~] = fileparts(outputFile);
    if ~isempty(outPath) && ~exist(outPath, 'dir')
        mkdir(outPath);
    end
    
    set(fig, 'Toolbar', 'none');
    saveas(fig, outputFile);
    fprintf('Figure saved successfully.\n');
    close(fig);
end

function res = computeResiduals(baseImageFile, regressionTarget, targetFile, dataDir)
    [~, baseName, ~] = fileparts(baseImageFile);
    
    if isempty(dataDir)
        [dataDir, ~, ~] = fileparts(which(mfilename));
        if isempty(dataDir)
            dataDir = pwd;
        end
    end
    
    templateCacheFile = fullfile(dataDir, ['template' upper(baseName(1)) lower(baseName(2:end)) '.mat']);
    
    if ~isfile(templateCacheFile)
        error('Template file not found: %s', templateCacheFile);
    end
    
    load(templateCacheFile);
    res = res;
end

function princ = computePrinc(regressionTarget, targetFile, baseName, dataDir)
    if isempty(dataDir)
        [dataDir, ~, ~] = fileparts(which(mfilename));
        if isempty(dataDir)
            dataDir = pwd;
        end
    end
    
    switch lower(regressionTarget)
        case 'morphology'
            morphFile = fullfile(dataDir, 'Lesions_PN_morphology.mat');
            if ~isfile(morphFile)
                morphFile = fullfile(dataDir, 'morphology.mat');
            end
            if ~isfile(morphFile)
                error('Morphology file not found: %s', morphFile);
            end
            load(morphFile);
            if ~exist('princ', 'var')
                error('princ variable not found in morphology file');
            end
            princ = princ(:,1);
            
        case 'targetfile'
            [~, targetName, ~] = fileparts(targetFile);
            targetCacheFile = fullfile(dataDir, ['template' upper(targetName(1)) lower(targetName(2:end)) '.mat']);
            if ~isfile(targetCacheFile)
                error('Target file not found: %s', targetCacheFile);
            end
            targetData = load(targetCacheFile);
            if isfield(targetData, 'res')
                princ = targetData.res(:,1);
            elseif isfield(targetData, 'princ')
                princ = targetData.princ;
            else
                error('Target file does not contain res or princ variable');
            end
            
        case 'baseimage'
            princ = [];
    end
end

function h = drawHalfViolin(xPos, data, side, faceColor, edgeColor, faceAlpha, maxWidth)
    % drawHalfViolin - Draw one half of a split violin plot
    %
    % Arguments:
    %   xPos      - x-axis center position for the violin
    %   data      - vector of values to estimate density from
    %   side      - 'left' or 'right'
    %   faceColor - RGB triplet for fill color
    %   edgeColor - RGB triplet for edge color
    %   faceAlpha - face transparency (0-1)
    %   maxWidth  - maximum half-width of the violin in x-axis units

    [f, xi] = ksdensity(data, 'NumPoints', 200);
    % Normalize density so the peak equals maxWidth
    f = f / max(f) * maxWidth;

    if strcmp(side, 'left')
        xCoords = [xPos - f, xPos * ones(1, length(xi))];
    else
        xCoords = [xPos * ones(1, length(xi)), xPos + fliplr(f)];
    end
    yCoords = [xi, fliplr(xi)];

    h = fill(xCoords, yCoords, faceColor, ...
        'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, 'LineWidth', 1);

    % Draw median line
    med = median(data);
    if strcmp(side, 'left')
        medWidth = interp1(xi, f, med, 'linear', 0);
        line([xPos - medWidth, xPos], [med, med], ...
            'Color', edgeColor, 'LineWidth', 2);
    else
        medWidth = interp1(xi, f, med, 'linear', 0);
        line([xPos, xPos + medWidth], [med, med], ...
            'Color', edgeColor, 'LineWidth', 2);
    end
end
