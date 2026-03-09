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
% This script contains the function for computing morphological         %
% features of lesions using SPM.                                        %
%                                                                       %
% input:                                                                %

% -lesionDir:  containing lesion files (default: '/path/to/Lesions')    %
% -dataDir: directory for masks and output (default: '/path/to/Data')   %
% -spmPath: Path to SPM installation (default: '/path/to/spm')          %
%                                                                       %
% last update: 2026/03/07.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preprocessingMorphology(varargin)
    % Parse command line arguments
    p = inputParser;
    addParameter(p, 'lesionDir', '/path/to/Lesions', @ischar);
    addParameter(p, 'dataDir', '/path/to/Data', @ischar);
    addParameter(p, 'spmPath', '/path/to/spm', @ischar);
    addParameter(p, 'help', false, @islogical);
    parse(p, varargin{:});
    
    if p.Results.help
        help preprocessingMorphology;
        return;
    end
    
    lesionDir = p.Results.lesionDir;
    dataDir = p.Results.dataDir;
    spmPath = p.Results.spmPath;
    
    fprintf('========================================\n');
    fprintf('PREPROCESSING MORPHOLOGY\n');
    fprintf('========================================\n');
    fprintf('Lesion directory: %s\n', lesionDir);
    fprintf('Data directory: %s\n', dataDir);
    fprintf('========================================\n\n');
    
    % Ensure SPM is in the MATLAB path
    if exist('spm', 'file') ~= 2
        fprintf('Adding SPM to path: %s\n', spmPath);
        addpath(spmPath);
    end
    
    % Define output file
    [~, lesionDirBase] = fileparts(lesionDir);
    morphFile = fullfile(dataDir, [lesionDirBase '_morphology.mat']);
    
    % Check if already exists
    if isfile(morphFile)
        fprintf('Morphology file already exists: %s\n', morphFile);
        load(morphFile);
        fprintf('Loaded existing morphology data\n');
        return;
    end
    
    % Get lesion files
    src = dir(fullfile(lesionDir, '*.nii'));
    if isempty(src)
        error('No .nii files found in lesion directory: %s', lesionDir);
    end
    fprintf('Found %d lesion files\n', length(src));
    
    fprintf('Computing morphology features...\n');
    centr = zeros(length(src), 3);
    princ = zeros(length(src), 3);
    orients = zeros(length(src), 3);
    surfa = zeros(length(src), 1);
    eigenvals = zeros(length(src), 3);
    
    for i = 1:length(src)
        tar = spm_read_vols(spm_vol(fullfile(src(i).folder, src(i).name)));
        stats = regionprops3(tar, 'Centroid', 'SurfaceArea', ...
            'PrincipalAxisLength', 'Orientation', 'Eigenvalues');
        
        centr(i,:) = stats.Centroid;
        princ(i,:) = stats.PrincipalAxisLength;
        orients(i,:) = stats.Orientation;
        surfa(i) = stats.SurfaceArea;
        eigenvals(i,:) = stats.EigenValues{1};
        
        if mod(i, 100) == 0
            fprintf('  Processed %d/%d lesions\n', i, length(src));
        end
    end
    
    % Save results (using original variable name 'princ' for compatibility)
    save(morphFile, 'princ', 'centr', 'orients', 'surfa', 'eigenvals', 'lesionDir');
    fprintf('\nSaved morphology file: %s\n', morphFile);
    fprintf('  - princ: Principal axis lengths (3 columns)\n');
    fprintf('  - centr: Centroid positions (3 columns)\n');
    fprintf('  - orients: Orientations (3 columns)\n');
    fprintf('  - surfa: Surface areas\n');
    fprintf('  - eigenvals: Eigenvalues (3 columns)\n');
    fprintf('\nPreprocessing complete!\n');
end

