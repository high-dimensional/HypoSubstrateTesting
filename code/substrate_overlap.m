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
% This script contains the function for computing the overlap of every  %
% lesion with a substrate mask using SPM image calculations. Volumes    %
% are given in 1000s of voxels: vIn (lesion inside the substrate), vOut %
% (lesion outside) and vTot (total lesion volume).                      %
%                                                                       %
% input:                                                                %
% -substrate: substrate name (<name>.nii in maskDir)                    %
% -maskDir: directory with the substrate masks                          %
% -lesionDir: directory with the lesion .nii files                      %
% -outDir: output directory                                             %
% -spmPath: path to the SPM installation                                %
%                                                                       %
%                                                                       %
% Author: Patrik Bey, patrik.bey@ucl.ac.uk                              %
%                                                                       %
%                                                                       %
% last update: 2026/11/10.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function substrate_overlap(varargin)
    p = inputParser;
    addParameter(p, 'substrate', '', @ischar);
    addParameter(p, 'maskDir', '', @ischar);
    addParameter(p, 'lesionDir', '', @ischar);
    addParameter(p, 'outDir', '', @ischar);
    addParameter(p, 'spmPath', '', @ischar);
    parse(p, varargin{:});
    substrate = p.Results.substrate;
    outDir = p.Results.outDir;
    if ~isempty(p.Results.spmPath) && exist('spm_imcalc', 'file') ~= 2
        addpath(p.Results.spmPath);
    end
    if ~isfolder(outDir), mkdir(outDir); end

    maskFile = fullfile(p.Results.maskDir, [substrate '.nii']);
    if ~isfile(maskFile), error('Substrate mask not found: %s', maskFile); end
    maskTotal = sum(spm_read_vols(spm_vol(maskFile)) > 0, 'all');

    src = dir(fullfile(p.Results.lesionDir, '*.nii'));
    if isempty(src), error('No .nii lesion files in %s', p.Results.lesionDir); end
    n = numel(src);
    lesionNames = {src.name}';
    vIn = nan(n, 1);
    vOut = nan(n, 1);

    % private scratch directory so that substrates can be processed in parallel
    scratch = tempname;
    mkdir(scratch);
    cleanup = onCleanup(@() rmdir(scratch, 's'));
    diffFile = fullfile(scratch, 'diff.nii');
    flags = struct('dmtx', 0, 'mask', 0, 'interp', 0, 'dtype', 4);

    for i = 1:n
        inputs = spm_vol(char(maskFile, fullfile(src(i).folder, src(i).name)));
        % +1: substrate voxel spared by the lesion, -1: lesion voxel outside the substrate
        spm_imcalc(inputs, diffFile, '(i1>0)-(i2>0)', flags);
        d = spm_read_vols(spm_vol(diffFile));
        vIn(i) = maskTotal - sum(d(:) > 0);
        vOut(i) = sum(d(:) < 0);
        if mod(i, 250) == 0, fprintf('  %s: %d/%d lesions\n', substrate, i, n); end
    end

    vTot = (vIn + vOut) / 1000;
    vIn = vIn / 1000;
    vOut = vOut / 1000;

    save(fullfile(outDir, ['overlap_' substrate '.mat']), 'lesionNames', 'vIn', 'vOut', 'vTot');
    fprintf('Saved overlap features for %s (%d lesions)\n', substrate, n);
end
