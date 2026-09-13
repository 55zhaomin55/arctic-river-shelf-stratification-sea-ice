function cfg = project_config()
%PROJECT_CONFIG Return repository-relative input and output paths.
%
% The repository is assumed to have the structure
%   code/       MATLAB source code
%   data/       processed inputs used by the figure scripts
%   output/     generated figures and tables
%
% No machine-specific absolute paths are used.

codeDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(codeDir);

cfg.rootDir   = rootDir;
cfg.codeDir   = codeDir;
cfg.figureDir = fullfile(codeDir,'figures');
cfg.funcDir   = fullfile(codeDir,'functions');
cfg.dataDir   = fullfile(rootDir,'data');
cfg.geoDir    = fullfile(cfg.dataDir,'geospatial');
cfg.outputDir = fullfile(rootDir,'output');

if ~exist(cfg.outputDir,'dir')
    mkdir(cfg.outputDir);
end
end
