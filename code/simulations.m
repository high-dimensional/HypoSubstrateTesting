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
% This script contains the function for running Monte Carlo simulations %
% to test whether substrate map associations with cognitive outcomes    %
% are spurious (driven by lesion volume confounds) or true.             %
%                                                                       %
% input:                                                                %
% -baseImageFile: Name of substrate image file (default: 'MDN.nii')     %
% -regressionTarget: Target type (default: 'morphology')                %
% -targetFile: Target file for overlap-based cognition                  %
% -lesionDir:  containing lesion files (default: '/path/to/Lesions')    %
% -dataDir: directory for masks and output (default: '/path/to/Data')   %
% -spmPath: Path to SPM installation (default: '/path/to/spm')          %
% -ncas: Number of cases per iteration (default: 300)                   %
% -nlow: Low noise level (default: 0.3)                                 %
% -nhigh: High noise level (default: 0.6)                               %
% -nits: Number of iterations (default: 5000)                           %
% -splitPercentile: Percentile threshold for balanced bootstrapping     %
% -outputDir: Directory for saving results                              %
% -figureFormat: Figure format (default: 'pdf')                         %
% -binSize: Noise-level binning (default: 0.05)                         %    
% -n_subset: Randomly subset the data to this N                         %
% -seed: Random seed for reproducibility (default: 42)                  %
%                                                                       %
% last update: 2026/03/07.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simulations(varargin)
    % Parse command line arguments
    p = inputParser;
    addParameter(p, 'baseImageFile', 'MDN.nii', @ischar);
    addParameter(p, 'regressionTarget', 'morphology', @ischar);
    addParameter(p, 'targetFile', 'RF.nii', @ischar);
    addParameter(p, 'lesionDir', '/path/to/Lesions', @ischar);
    addParameter(p, 'dataDir', '/path/to/Data', @ischar);
    addParameter(p, 'outputDir', '', @ischar);
    addParameter(p, 'spmPath', '/path/to/spm', @ischar);
    addParameter(p, 'ncas', 300, @isnumeric);
    addParameter(p, 'nlow', 0.3, @isnumeric);
    addParameter(p, 'nhigh', 0.6, @isnumeric);
    addParameter(p, 'nits', 5000, @isnumeric);
    addParameter(p, 'splitPercentile', 30, @isnumeric);
    addParameter(p, 'skipPreprocessing', false, @islogical);
    addParameter(p, 'figureFormat', 'pdf', @ischar);
    addParameter(p, 'binSize', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'n_subset', NaN, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'seed', 42, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'help', false, @islogical);
    parse(p, varargin{:});
    
    if p.Results.help
        help simulations;
        return;
    end
    
    % Extract parameters
    baseImageFile = p.Results.baseImageFile;
    regressionTarget = p.Results.regressionTarget;
    targetFile = p.Results.targetFile;
    lesionDir = p.Results.lesionDir;
    dataDir = p.Results.dataDir;
    outputDir = p.Results.outputDir;
    if isempty(outputDir)
        outputDir = dataDir;
    end
    spmPath = p.Results.spmPath;
    ncas = p.Results.ncas;
    nlow = p.Results.nlow;
    nhigh = p.Results.nhigh;
    nits = p.Results.nits;
    splitPercentile = p.Results.splitPercentile;
    skipPreprocessing = p.Results.skipPreprocessing;
    figureFormat = p.Results.figureFormat;
    binSize = p.Results.binSize;
    n_subset = p.Results.n_subset;
    seed = p.Results.seed;
    
    % Validate target type
    validTargets = {'morphology', 'targetfile'};
    if ~ismember(lower(regressionTarget), validTargets)
        error('Invalid regressionTarget. Must be one of: morphology, targetfile');
    end
    
    fprintf('========================================\n');
    fprintf('SIMULATION EXPERIMENTS\n');
    fprintf('========================================\n');
    fprintf('Base image (network mask): %s\n', baseImageFile);
    fprintf('Regression target: %s\n', regressionTarget);
    if strcmpi(regressionTarget, 'targetfile')
        fprintf('Target file: %s\n', targetFile);
    end
    fprintf('Lesion directory: %s\n', lesionDir);
    fprintf('Data directory: %s\n', dataDir);
    fprintf('Output directory: %s\n', outputDir);
    fprintf('Cases per iteration: %d\n', ncas);
    if isnan(n_subset)
        fprintf('Data subset: all (no subsetting)\n');
    else
        fprintf('Data subset: %d rows\n', n_subset);
    end
    fprintf('Split percentile: %d%% (low <= %dth pctl, high > %dth pctl)\n', splitPercentile, splitPercentile, splitPercentile);
    fprintf('Noise range: %.2f - %.2f\n', nlow, nhigh);
    fprintf('Iterations: %d\n', nits);
    fprintf('Random seed: %d\n', seed);
    fprintf('========================================\n\n');
    
    % Ensure SPM is in the MATLAB path
    if exist('spm', 'file') ~= 2
        fprintf('Adding SPM to path: %s\n', spmPath);
        addpath(spmPath);
    end
    
    % Generate output filenames
    [~, baseName, ~] = fileparts(baseImageFile);
    [~, targetName, ~] = fileparts(targetFile);
    
    % Consistent naming: base-{baseimagename}_{target}
    switch lower(regressionTarget)
        case 'targetfile'
            resultFile = fullfile(outputDir, sprintf('base-%s_target-%s.mat', baseName, targetName));
        case 'morphology'
            resultFile = fullfile(outputDir, sprintf('base-%s_morphology.mat', baseName));
    end
    
    % === LOAD/COMPUTE BASE IMAGE INTERSECTIONS =================================
    [~, baseName, ~] = fileparts(baseImageFile);
    templateCacheFile = fullfile(dataDir, ['template' upper(baseName(1)) lower(baseName(2:end)) '.mat']);
    
    if isfile(templateCacheFile)
        fprintf('Loading existing template file: %s\n', templateCacheFile);
        load(templateCacheFile);
    else
        if skipPreprocessing
            error('Template file not found and skipPreprocessing is true: %s', templateCacheFile);
        end
        fprintf('Template file not found. Running preprocessingFeatures...\n');
        preprocessingFeatures('baseImageFile', baseImageFile, 'lesionDir', lesionDir, ...
            'dataDir', dataDir, 'spmPath', spmPath);
        load(templateCacheFile);
    end
    
    % CRITICAL FIX: Save base image res before loading target data
    % This ensures res always contains base image overlaps for the simulations
    res_base = res;
    fprintf('Base image data loaded and preserved (res_base: %dx%d)\n', size(res_base, 1), size(res_base, 2));
    
    % === LOAD/COMPUTE TARGET-SPECIFIC DATA =====================================
    % Note: Target files also contain a variable called 'res', which would
    % overwrite our base image res if we don't preserve it first
    switch lower(regressionTarget)
        case 'morphology'
            % Use lesionDir basename for dynamic morphology filename
            [~, lesionDirBase] = fileparts(lesionDir);
            morphFile = fullfile(dataDir, [lesionDirBase '_morphology.mat']);
            % Fallback to legacy filename for backward compatibility
            if ~isfile(morphFile)
                legacyMorphFile = fullfile(dataDir, 'morphology.mat');
                if isfile(legacyMorphFile)
                    morphFile = legacyMorphFile;
                end
            end
            if isfile(morphFile)
                fprintf('Loading existing morphology file: %s\n', morphFile);
                load(morphFile);
            else
                if skipPreprocessing
                    error('Morphology file not found and skipPreprocessing is true: %s', morphFile);
                end
                fprintf('Morphology file not found. Running preprocessingMorphology...\n');
                preprocessingMorphology('lesionDir', lesionDir, 'dataDir', dataDir, 'spmPath', spmPath);
                load(morphFile);
            end
            % princ is already loaded from morphFile (contains principal axis lengths)
            if ~exist('princ', 'var')
                error('princ variable not found in morphology file');
            end
            % Use only the first column (primary principal axis length)
            princ = princ(:,1);
            
        case 'targetfile'
            [~, targetName, ~] = fileparts(targetFile);
            targetCacheFile = fullfile(dataDir, ['template' upper(targetName(1)) lower(targetName(2:end)) '.mat']);
            
            if isfile(targetCacheFile)
                fprintf('Loading existing target file: %s\n', targetCacheFile);
                % Load into temporary structure to avoid overwriting res_base
                targetData = load(targetCacheFile);
                % Check if res exists, otherwise use princ (backward compatibility)
                if isfield(targetData, 'res')
                    fprintf('Target file loaded (targetData.res: %dx%d)\n', size(targetData.res, 1), size(targetData.res, 2));
                    princ = targetData.res(:,1);
                elseif isfield(targetData, 'princ')
                    fprintf('Target file loaded (using princ field: %dx%d)\n', size(targetData.princ, 1), size(targetData.princ, 2));
                    princ = targetData.princ;
                else
                    error('Target file does not contain res or princ variable: %s', targetCacheFile);
                end
            else
                error('Target cache file not found: %s. Run preprocessingFeatures first.', targetCacheFile);
            end
            fprintf('princ extracted from target data (size: %dx%d)\n', size(princ, 1), size(princ, 2));
    end
    
    % Ensure res contains the base image data (not target data)
    res = res_base;
    clear res_base;  % Clean up temporary variable
    fprintf('res restored to base image data for simulations\n');
    
    % === SET RANDOM SEED (before any random operations) =========================
    rng(seed);
    fprintf('Random seed set to: %d\n', seed);

    % === OPTIONAL DATA SUBSETTING ==============================================
    if ~isnan(n_subset)
        nTotal = size(res, 1);
        if n_subset >= nTotal
            fprintf('n_subset (%d) >= total rows (%d). Using all data.\n', n_subset, nTotal);
        else
            fprintf('Subsetting data: selecting %d of %d rows at random...\n', n_subset, nTotal);
            subIdx = randperm(nTotal, n_subset)';
            res = res(subIdx, :);
            princ = princ(subIdx, :);
            fprintf('Data subsetted to %d rows.\n', n_subset);
        end
    end
    
    % === RUN SIMULATIONS =======================================================
    fprintf('\nRunning simulations...\n');
    nlev = linspace(nlow, nhigh, nits);
    result = zeros(2, nits);
    
    fprintf('Running %d iterations...\n', nits);
    
    for j = 1:length(nlev)
        if mod(j, 500) == 0
            fprintf('  Completed %d/%d iterations\n', j, nits);
        end
        
        [result(:, j), iterData] = runIteration(res, princ, regressionTarget, nlev(j), ncas, splitPercentile);
        
        % Capture scatter plot data from the first and last iterations
        if j == 1
            scatterData_low = iterData;
            scatterData_low.noiseLevel = nlev(j);
        end
        if j == length(nlev)
            scatterData_high = iterData;
            scatterData_high.noiseLevel = nlev(j);
        end
    end
    
    fprintf('Simulations complete.\n');
    
    % === SAVE RESULTS ==========================================================
    fprintf('\nSaving results...\n');
    % Include scatterData structs so plot_sim_results can plot actual
    % simulation iterations without regenerating noise or rerunning regressions
    save(resultFile, 'result', 'nlev', 'regressionTarget', 'baseImageFile', ...
        'targetFile', 'lesionDir', 'ncas', 'nits', 'splitPercentile', ...
        'n_subset', 'binSize', 'seed', ...
        'scatterData_low', 'scatterData_high');
    fprintf('Results saved to: %s\n', resultFile);
    
    % === GENERATE FIGURES (via plot_sim_results) ================================
    fprintf('\nGenerating figures...\n');
    [~, figureBaseName, ~] = fileparts(resultFile);
    figureOutputFile = fullfile(outputDir, [figureBaseName '_figures.' figureFormat]);
    
    try
        plot_sim_results(resultFile, 'outputFile', figureOutputFile, 'binSize', binSize);
        fprintf('Figures saved to: %s\n', figureOutputFile);
    catch ME
        fprintf('WARNING: Figure generation failed: %s\n', ME.message);
        for k = 1:length(ME.stack)
            fprintf('  In %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
    end
    
    % === EXPORT BINNED SIGNIFICANCE RATIOS CSV ==================================
    fprintf('\nComputing binned significance ratios...\n');
    csvFile = fullfile(outputDir, [figureBaseName '_significance.csv']);
    
    significanceThreshold = -log(0.05);
    [~, baseImageName, ~] = fileparts(baseImageFile);
    
    % Build bin edges from nlow to nhigh in steps of binSize
    binEdges = round(nlow:binSize:nhigh, 10);
    if binEdges(end) < nhigh
        binEdges(end+1) = round(binEdges(end) + binSize, 10);
    end
    nBins = length(binEdges) - 1;
    
    % Assign each iteration to a bin
    binIdx = discretize(nlev, binEdges, 'IncludedEdge', 'right');
    
    % Write CSV
    fid = fopen(csvFile, 'w');
    if fid == -1
        error('Could not open CSV file for writing: %s', csvFile);
    end
    fprintf(fid, 'filename,baseImage,regressionTarget,targetFile,noiseBinStart,noiseBinEnd,nIterations,pctSurviving_slopeComparison,pctSurviving_withVolumeCovariate\n');
    
    [~, resultFileName, resultFileExt] = fileparts(resultFile);
    
    for b = 1:nBins
        inBin = (binIdx == b);
        nInBin = sum(inBin);
        
        if nInBin == 0
            continue;
        end
        
        pctSurviving1 = 100 * sum(result(1, inBin) > significanceThreshold) / nInBin;
        pctSurviving2 = 100 * sum(result(2, inBin) > significanceThreshold) / nInBin;
        
        fprintf(fid, '"%s","%s","%s","%s",%.4f,%.4f,%d,%.2f,%.2f\n', ...
            [resultFileName resultFileExt], ...
            baseImageName, ...
            regressionTarget, ...
            targetFile, ...
            binEdges(b), ...
            binEdges(b+1), ...
            nInBin, ...
            pctSurviving1, ...
            pctSurviving2);
    end
    
    fclose(fid);
    fprintf('Significance ratios saved to: %s (%d bins)\n', csvFile, nBins);
    
    fprintf('\nSimulation complete!\n');
    fprintf('Target: %s\n', regressionTarget);
    fprintf('Base: %s\n', baseImageFile);
end


function [result, dataOut] = runIteration(res, princ, regressionTarget, noiseLevel, ncas, splitPercentile)
    % runIteration - Run a single simulation iteration with balanced subsampling
    %
    % Workflow:
    %   1. Compute cognitive scores for ALL samples
    %   2. Add noise to ALL cognitive scores
    %   3. Subsample to get balanced classes based on noisy cognitive scores:
    %      - "Low" group: samples with cog <= splitPercentile threshold
    %      - "High" group: samples with cog > splitPercentile threshold
    %      - Sample ncas/2 from each group
    %   4. Run regression analysis on subsampled data
    
    % Get total number of cases from input matrix
    nCases = size(res, 1);
    
    % Extract lesion overlap data
    BaseVol = res(:,1);    % Base image overlap volume
    nBaseVol = res(:,2);   % Non-base image overlap volume
    totvol = res(:,3);     % Total lesion volume

    % === STEP 1 & 2: Compute cognitive scores for ALL samples and add noise ===
    
    % Generate noise for ALL samples
    pd = makedist('Normal', 'mu', 0, 'sigma', noiseLevel);
    dis = random(pd, [nCases 1]);

    % Compute cognitive target for ALL samples
    switch lower(regressionTarget)
        case 'morphology'
            cog = 1 - rescale(princ(:));
        case 'targetfile'
            princ_temp = rescale(princ(:));  % Ensure column vector
            cog = 1 - princ_temp;
    end
    cog = cog(:);  % Ensure column vector
    
    % Add noise to ALL cognitive scores
    cog = cog + dis;

    % === STEP 3: Balanced subsampling based on noisy cognitive scores ===
    
    % Calculate percentile threshold
    threshold = prctile(cog, splitPercentile);
    
    % Identify low and high groups
    lowIdx = find(cog <= threshold);   % Indices of "low" cognitive group
    highIdx = find(cog > threshold);   % Indices of "high" cognitive group
    
    % Number of samples per group
    nPerGroup = floor(ncas / 2);
    
    % Handle edge cases where groups have fewer samples than requested
    if length(lowIdx) < nPerGroup
        warning('runIteration:InsufficientLowSamples', ...
            'Low group has only %d samples (requested %d). Using all available.', ...
            length(lowIdx), nPerGroup);
        nLowSample = length(lowIdx);
    else
        nLowSample = nPerGroup;
    end
    
    if length(highIdx) < nPerGroup
        warning('runIteration:InsufficientHighSamples', ...
            'High group has only %d samples (requested %d). Using all available.', ...
            length(highIdx), nPerGroup);
        nHighSample = length(highIdx);
    else
        nHighSample = nPerGroup;
    end
    
    % Randomly sample from each group
    lowSample = lowIdx(randperm(length(lowIdx), nLowSample));
    highSample = highIdx(randperm(length(highIdx), nHighSample));
    
    % Combine selected indices
    selectedIdx = [lowSample; highSample];
    
    % Create selection vector (for compatibility with existing visualization code)
    selvec = zeros(nCases, 1);
    selvec(selectedIdx) = 1;
    indvec = logical([selvec; selvec]);  % Duplicated for BaseVol and nBaseVol rows

    % === STEP 4: Build data structure and run regression ===
    
    % Prepare comparison data (duplicate rows for BaseVol vs nBaseVol comparison)
    ovrl = [BaseVol; nBaseVol];           % Overlap values
    voll = [totvol; totvol];              % Volume values (same for both)
    cogg = [cog; cog];                    % Cognitive scores (same for both)
    comp = zeros(size(cogg, 1), 1);       % Network indicator
    comp(1:nCases) = 1;                   % First half = base network (1), second half = non-base (0)

    % Create data table with selected samples only
    data = table(cogg, ovrl, voll, comp, 'VariableNames', {'Cognition', 'Overlap', 'Volume', 'Network'});
    data(~indvec, :) = [];

    % Fit model without volume covariate
    fit = fitlm(data, 'interactions', 'ResponseVar', 'Cognition', ...
        'PredictorVars', {'Overlap', 'Network'}, ...
        'CategoricalVar', {'Network'}, 'RobustOpts', 'off');

    anv = anova(fit);
    tmp = anv.pValue;
    result(1) = -log(tmp(3));

    % Fit model with volume covariate
    fit2 = fitlm(data, 'interactions', 'ResponseVar', 'Cognition', ...
        'PredictorVars', {'Overlap', 'Volume', 'Network'}, ...
        'CategoricalVar', {'Network'}, 'RobustOpts', 'off');

    anv2 = anova(fit2);
    tmp = anv2.pValue;
    result(2) = -log(tmp(5));

    % Output data structure (for visualization)
    % Note: LinearModel objects don't serialize reliably across MATLAB
    % versions, so we save the filtered data table instead and refit at
    % plot time (cheap for a single iteration).
    dataOut.ovrl = ovrl;
    dataOut.cogg = cogg;
    dataOut.comp = comp;
    dataOut.indvec = indvec;
    dataOut.voll = voll;
    dataOut.data = data;
end
