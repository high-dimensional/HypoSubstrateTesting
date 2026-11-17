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
% This script contains the function for the regionally weighted lesion  %
% volume-deficit regression (RLDR) of one simulated study. Each patient %
% contributes a target row (Overlap = vIn, Compartment = 1) and a non-  %
% target row (Overlap = vOut, Compartment = 0) with identical behaviour %
% and volume. Model 1: Behaviour ~ Overlap * Compartment. Model 2: all  %
% pairwise terms of Overlap, Compartment and Volume.                    %
%                                                                       %
% input:                                                                %
% -behaviour, -vIn, -vOut, -vTot: column vectors over t patients        %
% output: p (two-sided, Type III), beta and tstat of Overlap:Compartment%
% term for [Model 1, Model 2], nlogp = -ln(p), and the long-format table%
%                                                                       %
%                                                                       %
% Author: Patrik Bey, patrik.bey@ucl.ac.uk                              %
%                                                                       %
%                                                                       %
% last update: 2026/11/10.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fit_rldr(behaviour, vIn, vOut, vTot)
    n = numel(behaviour);
    tbl = table([vIn; vOut], categorical([ones(n, 1); zeros(n, 1)]), [vTot; vTot], ...
        [behaviour; behaviour], 'VariableNames', {'Overlap', 'Compartment', 'Volume', 'Behaviour'});

    m1 = fitlm(tbl, 'interactions', 'ResponseVar', 'Behaviour', ...
        'PredictorVars', {'Overlap', 'Compartment'}, 'CategoricalVars', 'Compartment');
    m2 = fitlm(tbl, 'interactions', 'ResponseVar', 'Behaviour', ...
        'PredictorVars', {'Overlap', 'Compartment', 'Volume'}, 'CategoricalVars', 'Compartment');

    [p1, b1, t1] = interactionStats(m1);
    [p2, b2, t2] = interactionStats(m2);
    out.p = [p1, p2];
    out.beta = [b1, b2];
    out.tstat = [t1, t2];
    out.nlogp = -log(out.p(:));
    out.tbl = tbl;
end


function [p, beta, tstat] = interactionStats(mdl)
    a = anova(mdl, 'component', 3);   % Type III sums of squares
    row = strcmp(a.Properties.RowNames, 'Overlap:Compartment');
    if nnz(row) ~= 1
        error('Overlap:Compartment term not found in ANOVA table');
    end
    p = a.pValue(row);

    % reference level '0' (non-target): beta = slope(target) - slope(non-target)
    names = mdl.CoefficientNames;
    idx = find(startsWith(names, 'Overlap:Compartment') & ~contains(names, 'Volume'));
    if numel(idx) ~= 1
        error('Could not identify the Overlap:Compartment coefficient');
    end
    beta = mdl.Coefficients.Estimate(idx);
    tstat = mdl.Coefficients.tStat(idx);

    % single-degree-of-freedom term: Type III F-test and coefficient t-test coincide
    pCoef = mdl.Coefficients.pValue(idx);
    if abs(pCoef - p) > 1e-6 * p && abs(pCoef - p) > 1e-12
        error('Coefficient p (%g) and ANOVA p (%g) disagree', pCoef, p);
    end
end
