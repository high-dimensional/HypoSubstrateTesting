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
% This script contains the function for preprocessing features by       %
% computing lesion overlays with a substrate image using SPM.           %
%                                                                       %
% input:                                                                %
% -baseImageFile: Name of substrate image file (default: 'MDN.nii')     %
% -lesionDir:  containing lesion files (default: '/path/to/Lesions')    %
% -dataDir: directory for masks and output (default: '/path/to/Data')   %
% -spmPath: Path to SPM installation (default: '/path/to/spm')          %
%                                                                       %
% IMPORTANT                                                             %
% All lesion masks need to be binary and in the same space as the       %
% base image for example MNI152_2mm as in the current study.            %
% Files need to be in .nii format for use with SPM.                     %
%                                                                       %
% last update: 2026/03/07.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function preprocessingFeatures(varargin)
    % Parse command line arguments
    p = inputParser;
    addParameter(p, 'baseImageFile', 'MDN.nii', @ischar);
    addParameter(p, 'lesionDir', '/path/to/Lesions', @ischar);
    addParameter(p, 'dataDir', '/path/to/Data', @ischar);
    addParameter(p, 'spmPath', '/path/to/spm', @ischar);
    addParameter(p, 'help', false, @islogical);
    parse(p, varargin{:});
    
    if p.Results.help
        help preprocessingFeatures;
        return;
    end
    
    baseImageFile = p.Results.baseImageFile;
    lesionDir = p.Results.lesionDir;
    dataDir = p.Results.dataDir;
    spmPath = p.Results.spmPath;
    
    fprintf('========================================\n');
    fprintf('PREPROCESSING FEATURES\n');
    fprintf('========================================\n');
    fprintf('Base image: %s\n', baseImageFile);
    fprintf('Lesion directory: %s\n', lesionDir);
    fprintf('Data directory: %s\n', dataDir);
    fprintf('========================================\n\n');
    
    % Ensure SPM is in the MATLAB path
    if exist('spm', 'file') ~= 2
        fprintf('Adding SPM to path: %s\n', spmPath);
        addpath(spmPath);
    end
    
    % Load substrate image
    mskPath = fullfile(dataDir, 'Masks', baseImageFile);
    if ~isfile(mskPath)
        error('Base image file not found: %s', mskPath);
    end
    
    msk = spm_read_vols(spm_vol(mskPath));
    msktot = sum(msk(:) > 0);
    fprintf('Loaded base image: %d voxels in mask\n', msktot);
    
    % Get lesion files
    src = dir(fullfile(lesionDir, '*.nii'));
    if isempty(src)
        error('No .nii files found in lesion directory: %s', lesionDir);
    end
    fprintf('Found %d lesion files\n', length(src));
    
    % Create cache filename based on base image
    [~, baseName, ~] = fileparts(baseImageFile);
    templateCacheFile = fullfile(dataDir, ['template' upper(baseName(1)) lower(baseName(2:end)) '.mat']);
    
    fprintf('Computing intersections...\n');
    res = zeros(length(src), 3);
    
    for i = 1:length(src)
        matlabbatch{1}.spm.util.imcalc.input = {
            mskPath
            [fullfile(src(i).folder, src(i).name) ',1']
        };
        matlabbatch{1}.spm.util.imcalc.output = ['tmp-', baseImageFile];
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1-i2';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch);
        clear matlabbatch
        
        tar = spm_read_vols(spm_vol(['tmp-',baseImageFile]));
        res(i,1) = msktot - sum(tar(:) > 0);    % lesioned base image volume
        res(i,2) = sum(tar(:) < 0);              % lesion volume outside base image
        res(i,3) = res(i,2) + res(i,1);          % total lesion volume
        
        if mod(i, 100) == 0
            fprintf('  Processed %d/%d lesions\n', i, length(src));
        end
    end
    
    % Rescale to thousands of voxels
    res(:,1) = res(:,1) ./ 1000;  % lesioned base image volume
    res(:,2) = res(:,2) ./ 1000;  % lesion volume outside base image
    res(:,3) = res(:,3) ./ 1000;  % total lesion volume
    
    % Save results
    save(templateCacheFile, 'res', 'baseImageFile', 'lesionDir');
    fprintf('\nSaved template file: %s\n', templateCacheFile);
    fprintf('  - res(:,1): Lesioned base image volume (in 1000s voxels)\n');
    fprintf('  - res(:,2): Lesion volume outside base image (in 1000s voxels)\n');
    fprintf('  - res(:,3): Total lesion volume (in 1000s voxels)\n');
    fprintf('\nPreprocessing complete!\n');
end
