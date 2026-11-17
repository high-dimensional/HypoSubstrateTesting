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
% This script contains the function for computing the first principal   %
% axis length of each lesion (regionprops3), the lesion size feature of %
% the global scenario. Multi-component lesions are treated as one       %
% object. Can be run in parts (part/nParts) and merged (merge).         %
%                                                                       %
% input:                                                                %
% -lesionDir: directory with the lesion .nii files                      %
% -outDir: output directory                                             %
% -spmPath: path to the SPM installation                                %
% -part, -nParts: process lesion part of nParts                         %
% -merge: merge the parts into lesion_features.mat                      %
%                                                                       %
%                                                                       %
% Author: Patrik Bey, patrik.bey@ucl.ac.uk                              %
%                                                                       %
%                                                                       %
% last update: 2026/11/10.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lesion_features(varargin)
    p = inputParser;
    addParameter(p, 'lesionDir', '', @ischar);
    addParameter(p, 'outDir', '', @ischar);
    addParameter(p, 'spmPath', '', @ischar);
    addParameter(p, 'part', 1, @isnumeric);
    addParameter(p, 'nParts', 1, @isnumeric);
    addParameter(p, 'merge', false, @islogical);
    parse(p, varargin{:});
    lesionDir = p.Results.lesionDir;
    outDir = p.Results.outDir;
    part = p.Results.part;
    nParts = p.Results.nParts;
    if ~isempty(p.Results.spmPath) && exist('spm_vol', 'file') ~= 2
        addpath(p.Results.spmPath);
    end
    if ~isfolder(outDir), mkdir(outDir); end

    src = dir(fullfile(lesionDir, '*.nii'));
    if isempty(src), error('No .nii lesion files in %s', lesionDir); end
    n = numel(src);
    lesionNames = {src.name}';

    if p.Results.merge
        principalAxis = nan(n, 1);
        for k = 1:nParts
            P = load(fullfile(outDir, sprintf('lesion_features_part%d.mat', k)));
            assert(isequal(P.lesionNames, lesionNames), 'Lesion list changed between parts');
            principalAxis(P.idx) = P.axisPart;
        end
        assert(~any(isnan(principalAxis)), 'Missing principal axis values after merge');
        save(fullfile(outDir, 'lesion_features.mat'), 'lesionNames', 'principalAxis');
        fprintf('Merged %d parts: %d lesion features saved to %s\n', nParts, n, outDir);
        return;
    end

    idx = (part:nParts:n)';
    axisPart = nan(numel(idx), 1);
    for t = 1:numel(idx)
        i = idx(t);
        % numeric 0/1 input is treated as a label matrix: all components form one object
        vol = double(spm_read_vols(spm_vol(fullfile(src(i).folder, src(i).name))) > 0);
        stats = regionprops3(vol, 'PrincipalAxisLength');
        if height(stats) ~= 1
            error('Lesion %s: expected one region, got %d', src(i).name, height(stats));
        end
        axisPart(t) = stats.PrincipalAxisLength(1);
        if mod(t, 100) == 0, fprintf('  part %d: %d/%d lesions\n', part, t, numel(idx)); end
    end

    if nParts == 1
        principalAxis = axisPart;
        save(fullfile(outDir, 'lesion_features.mat'), 'lesionNames', 'principalAxis');
    else
        save(fullfile(outDir, sprintf('lesion_features_part%d.mat', part)), 'lesionNames', 'idx', 'axisPart');
    end
    fprintf('Part %d/%d done (%d lesions)\n', part, nParts, numel(idx));
end
