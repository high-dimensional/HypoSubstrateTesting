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
% This script contains the function for running the semi-synthetic      %
% Monte Carlo simulation. Behaviour is generated as y = 1 - rescale(x)  %
% + noise with noise ~ N(0, sigma^2), sigma rising linearly from nlow   %
% to nhigh, where x is the first principal axis length of the lesion    %
% (generator 'global') or its overlap with a substrate (generator =   %
% substrate name).                                                    %
% Per iteration all lesions are split at the splitPercentile of the     %
% noisy behaviour, nPerGroup lesions are drawn from each side and RLDR  %
% (fit_rldr) is fitted with the hypothesized compartment as target.     %
% Significant interactions (two-sided p < 0.05) are classified by the   %
% sign of beta: H-implicated (sign(beta) == expectedSign) or            %
% ~H-implicated. H-implicated p-values are reported as false positives  %
%                                                                       %
% input:                                                                %
% -hypothesized: compartment tested as target (default: 'MDN')          %
% -generator: 'global' or the substrate-based behaviour                 %
% -featureDir: directory with lesion_features.mat and overlap_<name>.mat%
% -outDir: output directory                                             %
% -nits: number of iterations (default: 5000)                           %
% -nlow, -nhigh: noise range (default: 0.3 - 0.6)                       %
% -nPerGroup: lesions drawn per side (default: 150)                     %
% -splitPercentile: percentile of the split (default: 30)               %
% -binWidth: noise bin width of the summary (default: 0.05)             %
% -seed: random seed (default: 42)                                      %
% -expectedSign: sign of beta that favours H                            %
% (default: -1, performance score)                                      %
%                                                                       %
%                                                                       %
% Author: Patrik Bey, patrik.bey@ucl.ac.uk                              %
%                                                                       %
%                                                                       %
% last update: 2026/11/10.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_rldr_simulation(varargin)
    p = inputParser;
    addParameter(p, 'hypothesized', 'MDN', @ischar);
    addParameter(p, 'generator', 'global', @ischar);
    addParameter(p, 'featureDir', '', @ischar);
    addParameter(p, 'outDir', '', @ischar);
    addParameter(p, 'nits', 5000, @isnumeric);
    addParameter(p, 'nlow', 0.3, @isnumeric);
    addParameter(p, 'nhigh', 0.6, @isnumeric);
    addParameter(p, 'nPerGroup', 150, @isnumeric);
    addParameter(p, 'splitPercentile', 30, @isnumeric);
    addParameter(p, 'binWidth', 0.05, @isnumeric);
    addParameter(p, 'seed', 42, @isnumeric);
    addParameter(p, 'expectedSign', -1, @(x) isscalar(x) && any(x == [-1 1]));
    parse(p, varargin{:});
    hyp = p.Results.hypothesized;
    gen = p.Results.generator;
    featureDir = p.Results.featureDir;
    outDir = p.Results.outDir;
    nits = p.Results.nits;
    nlow = p.Results.nlow;
    nhigh = p.Results.nhigh;
    nPerGroup = p.Results.nPerGroup;
    splitPercentile = p.Results.splitPercentile;
    binWidth = p.Results.binWidth;
    seed = p.Results.seed;
    expectedSign = p.Results.expectedSign;
    if ~isfolder(outDir), mkdir(outDir); end

    scenario = sprintf('hyp-%s_true-%s', hyp, gen);
    fprintf('Scenario %s\n', scenario);

    % ---- load features (all files list the lesions in the same order) ----
    F = load(fullfile(featureDir, 'lesion_features.mat'));
    H = load(fullfile(featureDir, ['overlap_' hyp '.mat']));
    assert(isequal(F.lesionNames, H.lesionNames), 'Lesion order differs between feature files');
    if strcmp(gen, 'global')
        x = F.principalAxis;
    else
        G = load(fullfile(featureDir, ['overlap_' gen '.mat']));
        assert(isequal(F.lesionNames, G.lesionNames), 'Lesion order differs between feature files');
        x = G.vIn;
    end
    nLesions = numel(x);
    xs = rescale(x(:));

    % ---- Monte Carlo loop ----
    rng(seed);
    sigma = nlow + (0:nits - 1) * (nhigh - nlow) / (nits - 1);
    nlogp = nan(2, nits);
    pval = nan(2, nits);
    beta = nan(2, nits);
    tstat = nan(2, nits);
    examples = struct('iteration', {}, 'sigma', {}, 'p', {}, 'beta', {}, 'tbl', {});

    for j = 1:nits
        y = 1 - xs + sigma(j) * randn(nLesions, 1);

        thr = prctile(y, splitPercentile);
        lowIdx = find(y <= thr);
        highIdx = find(y > thr);
        sel = [lowIdx(randperm(numel(lowIdx), nPerGroup)); ...
               highIdx(randperm(numel(highIdx), nPerGroup))];

        fit = fit_rldr(y(sel), H.vIn(sel), H.vOut(sel), H.vTot(sel));
        nlogp(:, j) = fit.nlogp;
        pval(:, j) = fit.p(:);
        beta(:, j) = fit.beta(:);
        tstat(:, j) = fit.tstat(:);

        if j == 1 || j == nits
            examples(end + 1) = struct('iteration', j, 'sigma', sigma(j), ...
                'p', fit.p, 'beta', fit.beta, 'tbl', fit.tbl); %#ok<AGROW>
        end
        if mod(j, 500) == 0, fprintf('  %d/%d iterations\n', j, nits); end
    end

    % ---- directional classification ----
    alpha = 0.05;
    isHImplicated = sign(beta) == expectedSign;
    pOne = pval / 2;
    pOne(~isHImplicated) = 1 - pval(~isHImplicated) / 2;
    sig = pval < alpha;

    % ---- save results ----
    params = struct('hypothesized', hyp, 'generator', gen, 'nits', nits, 'nlow', nlow, ...
        'nhigh', nhigh, 'nPerGroup', nPerGroup, 'splitPercentile', splitPercentile, ...
        'binWidth', binWidth, 'seed', seed, 'nLesions', nLesions, 'expectedSign', expectedSign);
    save(fullfile(outDir, [scenario '.mat']), 'nlogp', 'pval', 'beta', 'tstat', 'pOne', 'sigma', 'examples', 'params');

    iterations = table((1:nits)', sigma(:), nlogp(1, :)', nlogp(2, :)', pval(1, :)', pval(2, :)', ...
        beta(1, :)', beta(2, :)', tstat(1, :)', tstat(2, :)', pOne(1, :)', pOne(2, :)', ...
        'VariableNames', {'iteration', 'sigma', 'nlogp_unadjusted', 'nlogp_adjusted', ...
        'p_unadjusted', 'p_adjusted', 'beta_unadjusted', 'beta_adjusted', ...
        't_unadjusted', 't_adjusted', 'pOneSided_unadjusted', 'pOneSided_adjusted'});
    writetable(iterations, fullfile(outDir, [scenario '_iterations.csv']));

    % ---- percentage of iterations per noise bin ----
    edges = round(nlow:binWidth:nhigh, 10);
    bin = discretize(sigma, edges);
    nBins = numel(edges) - 1;
    pct = @(x, in) 100 * mean(x(in));
    rows = zeros(nBins, 11);
    for b = 1:nBins
        in = bin == b;
        for m = 1:2
            rows(b, 0 + m) = pct(sig(m, :), in);                         % any direction
            rows(b, 2 + m) = pct(sig(m, :) & isHImplicated(m, :), in);   % H-implicated
            rows(b, 4 + m) = pct(sig(m, :) & ~isHImplicated(m, :), in);  % ~H-implicated
            rows(b, 6 + m) = pct(~sig(m, :), in);                        % not significant
            rows(b, 8 + m) = pct(pOne(m, :) < alpha, in);                % one-sided
        end
        rows(b, 11) = nnz(in);
    end
    summary = table(repmat(string(scenario), nBins, 1), repmat(string(hyp), nBins, 1), ...
        repmat(string(gen), nBins, 1), edges(1:end - 1)', edges(2:end)', rows(:, 11), ...
        rows(:, 1), rows(:, 2), rows(:, 3), rows(:, 4), rows(:, 5), rows(:, 6), ...
        rows(:, 7), rows(:, 8), rows(:, 9), rows(:, 10), ...
        'VariableNames', {'scenario', 'hypothesized', 'trueGenerator', 'noiseBinStart', ...
        'noiseBinEnd', 'nIterations', 'pctSignificant_unadjusted', 'pctSignificant_adjusted', ...
        'pctTargetFavouring_unadjusted', 'pctTargetFavouring_adjusted', ...
        'pctNonTargetFavouring_unadjusted', 'pctNonTargetFavouring_adjusted', ...
        'pctNotSignificant_unadjusted', 'pctNotSignificant_adjusted', ...
        'pctOneSided_unadjusted', 'pctOneSided_adjusted'});
    writetable(summary, fullfile(outDir, [scenario '_summary.csv']));
    disp(summary);
end
