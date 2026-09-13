%% ========================================================================
% prepare_JGROceans_public_inputs_local.m
%
% LOCAL DATA PREPARATION TEMPLATE
%
% Purpose
% -------
% Read the author's working data files from the local computer and create
% the standardized processed input files required by the public JGR: Oceans
% reproducibility package.
%
% IMPORTANT
% ---------
% 1. Edit the source-directory paths below for your local machine.
% 2. The generated files in publication_inputs/data/ are the files intended
%    for the public reproducibility package.
% 3. No third-party raw ORAS5/ERA5/NSIDC/ArcticGRO products are duplicated
%    here unless they are already processed inputs actually used by the
%    figure scripts.
%
% Generated public filenames
% --------------------------
% annual_river_shelf_metrics_1979_2024.mat
% arctic_hsi_climatology_1979_2024.mat
% marginal_sea_salinity_profiles_1979_2024.mat
% runoff_salinity_distance_bins_1979_2024.mat
% domain_sensitivity_hsi_sic_1979_2024.mat
% sea_ice_monthly_1979_2024.mat
% sea_ice_grid_coordinates.mat
% wind_speed_JAS_1979_2024.mat
% surface_air_temperature_JAS_1979_2024.mat
% arctic_oscillation_monthly_1979_2024.mat
% discharge_<river>_1972_2024.mat (six files)
% geospatial/basin_<river>.mat (six files)
% geospatial/arctic_marginal_seas.*
%% ========================================================================

clear; clc;

%% ========================================================================
% 1. USER SETTINGS
%% ========================================================================

% ---- Main working directories on this computer ----
SALINITY_DIR = '<EDIT_ME>/salinity';
SEAICE_DIR   = '<EDIT_ME>/concentration';
ATMOS_DIR    = '<EDIT_ME>/seaice';
RIVER_DIR    = '<EDIT_ME>/river';
BASIN_DIR    = '<EDIT_ME>/basin_boundaries';

% ---- Output directory ----
% By default, create a clean data folder beside this script.
% You may instead replace OUTPUT_ROOT with the public repository root.
SCRIPT_DIR  = fileparts(mfilename('fullpath'));
OUTPUT_ROOT = fullfile(SCRIPT_DIR,'publication_inputs');
DATA_OUT    = fullfile(OUTPUT_ROOT,'data');
GEO_OUT     = fullfile(DATA_OUT,'geospatial');

% Overwrite existing generated files?
OVERWRITE = true;

% ---- Optional manual source paths ----
% Leave as '' to use the automatic search implemented below.
MANUAL_MAJOR_RIVER_SHP = '';

% If automatic discharge-file discovery fails for any river, provide the
% exact file path here. Otherwise leave ''.
MANUAL_DISCHARGE = struct( ...
    'mackenzie','', ...
    'yukon','', ...
    'kolyma','', ...
    'lena','', ...
    'yenisei','', ...
    'ob','');

%% ========================================================================
% 2. STANDARD SOURCE FILES ALREADY KNOWN FROM THE WORKING DIRECTORY
%% ========================================================================

src = struct();

src.annualMetrics = fullfile(SALINITY_DIR,'HSI.mat');
src.hsiClimatology = first_existing({ ...
    fullfile(SALINITY_DIR,'ORAS5_HSI_JAS_climatology_1979_2024.mat'), ...
    fullfile(SALINITY_DIR,'ORAS5_HSI_JAS_1979_2024.mat')});

src.salinityProfiles = first_existing({ ...
    fullfile(SALINITY_DIR,'ORAS5_marginal_sea_monthly_profiles_1979_2024.mat'), ...
    fullfile(SALINITY_DIR,'marginal_sea_monthly_profiles_1979_2024.mat')});

src.runoffSalinityBins = first_existing({ ...
    fullfile(SALINITY_DIR,'ORAS5_S15_distance_bins_1979_2024.mat'), ...
    fullfile(SALINITY_DIR,'runoff_salinity_distance_bins_1979_2024.mat')});

src.domainSensitivity = first_existing({ ...
    fullfile(SALINITY_DIR,'Figure6_HSI_SIC_domain_sensitivity_final.mat'), ...
    fullfile(SALINITY_DIR,'domain_sensitivity_hsi_sic_1979_2024.mat')});

src.seaIceMonthly = fullfile(SEAICE_DIR,'sic.mat');
src.seaIceCoords  = fullfile(SEAICE_DIR,'lon_lat.mat');
src.wind          = fullfile(ATMOS_DIR,'wind','wind.mat');
src.sat           = fullfile(ATMOS_DIR,'wind','SAT.mat');
src.ao            = fullfile(ATMOS_DIR,'AO_index.mat');

src.marginalSeaShp = fullfile(ATMOS_DIR,'World_Seas_IHO_v3','Arctic_Marginal_Seas.shp');

% Basin-boundary source files
src.basins = struct( ...
    'mackenzie', fullfile(BASIN_DIR,'maken','maken.mat'), ...
    'yukon',     fullfile(BASIN_DIR,'yukon','yukon.mat'), ...
    'kolyma',    fullfile(BASIN_DIR,'kolyma','kolyma.mat'), ...
    'lena',      fullfile(BASIN_DIR,'lena','lena.mat'), ...
    'yenisei',   fullfile(BASIN_DIR,'yenisei','yenisei.mat'), ...
    'ob',        fullfile(BASIN_DIR,'ob','ob.mat'));

%% ========================================================================
% 3. CREATE OUTPUT DIRECTORIES
%% ========================================================================

ensure_dir(OUTPUT_ROOT);
ensure_dir(DATA_OUT);
ensure_dir(GEO_OUT);

logFile = fullfile(OUTPUT_ROOT,'data_preparation_log.txt');
fid = fopen(logFile,'w');
if fid < 0
    error('Cannot create log file: %s',logFile);
end
cleanupLog = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'JGR: Oceans reproducibility input preparation log\n');
fprintf(fid,'Generated: %s\n\n',char(datetime('now')));

fprintf('\n============================================================\n');
fprintf('Preparing standardized publication input files\n');
fprintf('Output: %s\n',DATA_OUT);
fprintf('============================================================\n\n');

%% ========================================================================
% 4. FINAL ANNUAL RIVER-SHELF METRICS (Figures 8-11, S2-S3)
%% ========================================================================

dst = fullfile(DATA_OUT,'annual_river_shelf_metrics_1979_2024.mat');
D = load_required(src.annualMetrics,'Final HSI.mat');

OUT = struct();

keep = { ...
    'systems','years', ...
    'HSI_JAS','HSI_JAS_dt','HSI_JAS_z', ...
    'SIC_JAS','SIC_JAS_dt','SIC_JAS_z', ...
    'Q_MaySep_km3','Q_MaySep_dt','Q_MaySep_z', ...
    'R_HSI_Q','P_HSI_Q','R_HSI_SIC','P_HSI_SIC'};

OUT = copy_fields_if_present(D,OUT,keep);

% Add standardized system names if not already available.
if ~isfield(OUT,'systems')
    OUT.systems = { ...
        'Mackenzie-Beaufort', ...
        'Yukon-Chukchi', ...
        'Kolyma-East Siberian', ...
        'Lena-Laptev', ...
        'Yenisei-Kara', ...
        'Ob-Kara'};
end
if ~isfield(OUT,'years')
    OUT.years = (1979:2024)';
end

required = {'HSI_JAS','HSI_JAS_dt','SIC_JAS','SIC_JAS_dt'};
assert_fields(OUT,required,'annual_river_shelf_metrics_1979_2024.mat');
save_struct(dst,OUT,OVERWRITE);
record(fid,src.annualMetrics,dst);

%% ========================================================================
% 5. ARCTIC HSI CLIMATOLOGY (Figure 3a)
%% ========================================================================

dst = fullfile(DATA_OUT,'arctic_hsi_climatology_1979_2024.mat');
D = load_required(src.hsiClimatology,'HSI climatology');
OUT = struct();
OUT = copy_fields_if_present(D,OUT,{ ...
    'lonGrid','latGrid','HSI_JAS_clim','HSI_valid_years','years', ...
    'depthRead','depth','g','beta0','Hmax','minValidYears'});
assert_fields(OUT,{'lonGrid','latGrid','HSI_JAS_clim'}, ...
    'arctic_hsi_climatology_1979_2024.mat');
save_struct(dst,OUT,OVERWRITE);
record(fid,src.hsiClimatology,dst);

%% ========================================================================
% 6. MARGINAL-SEA MONTHLY SALINITY PROFILES (Figure 2)
%% ========================================================================

dst = fullfile(DATA_OUT,'marginal_sea_salinity_profiles_1979_2024.mat');
D = load_required(src.salinityProfiles,'Marginal-sea salinity profiles');
OUT = struct();
OUT = copy_fields_if_present(D,OUT,{ ...
    'SAL_MONTHLY','NUM_MONTHLY','DEN_MONTHLY','TIME_COMMON','depth', ...
    'sea_names','seaNames','years','months'});

% Standardize sea_names if source uses seaNames.
if ~isfield(OUT,'sea_names') && isfield(OUT,'seaNames')
    OUT.sea_names = OUT.seaNames;
end
assert_fields(OUT,{'SAL_MONTHLY','TIME_COMMON','depth'}, ...
    'marginal_sea_salinity_profiles_1979_2024.mat');
save_struct(dst,OUT,OVERWRITE);
record(fid,src.salinityProfiles,dst);

%% ========================================================================
% 7. RUNOFF-RELATED SALINITY DISTANCE BINS (Figure 5 and Figure S1)
%% ========================================================================

dst = fullfile(DATA_OUT,'runoff_salinity_distance_bins_1979_2024.mat');
D = load_required(src.runoffSalinityBins,'Runoff-salinity distance bins');
OUT = struct();
OUT = copy_fields_if_present(D,OUT,{ ...
    'SAL_DIST','N_VALID_CELLS','TIME_COMMON', ...
    'dist_edges','dist_centers','distance_edges','distance_centers', ...
    'lon','lat','years','months'});

% Standardize alternative distance names.
if ~isfield(OUT,'dist_edges') && isfield(OUT,'distance_edges')
    OUT.dist_edges = OUT.distance_edges;
end
if ~isfield(OUT,'dist_centers') && isfield(OUT,'distance_centers')
    OUT.dist_centers = OUT.distance_centers;
end
assert_fields(OUT,{'SAL_DIST','TIME_COMMON'}, ...
    'runoff_salinity_distance_bins_1979_2024.mat');
save_struct(dst,OUT,OVERWRITE);
record(fid,src.runoffSalinityBins,dst);

%% ========================================================================
% 8. DOMAIN-SENSITIVITY RESULTS (Figure 6)
%% ========================================================================

dst = fullfile(DATA_OUT,'domain_sensitivity_hsi_sic_1979_2024.mat');
D = load_required(src.domainSensitivity,'Figure 6 domain sensitivity');
OUT = struct();
OUT = copy_fields_if_present(D,OUT,{ ...
    'Rmat','Pmat','Nmat','TS_HSI','TS_SIC', ...
    'radii_km','center_offset','riverNames', ...
    'mouth_lon','mouth_lat','yrs_all','g','beta0','Hmax'});
assert_fields(OUT,{'Rmat','Pmat','Nmat','radii_km','center_offset'}, ...
    'domain_sensitivity_hsi_sic_1979_2024.mat');
save_struct(dst,OUT,OVERWRITE);
record(fid,src.domainSensitivity,dst);

%% ========================================================================
% 9. MONTHLY SEA ICE AND GRID COORDINATES (Figures 3 and 9)
%% ========================================================================

% Monthly SIC
dst = fullfile(DATA_OUT,'sea_ice_monthly_1979_2024.mat');
D = load_required(src.seaIceMonthly,'Monthly SIC');
if ~isfield(D,'sic')
    error('Variable "sic" not found in %s',src.seaIceMonthly);
end
OUT = struct('sic',D.sic);
save_struct(dst,OUT,OVERWRITE);
record(fid,src.seaIceMonthly,dst);

% SIC coordinates
dst = fullfile(DATA_OUT,'sea_ice_grid_coordinates.mat');
D = load_required(src.seaIceCoords,'SIC grid coordinates');
assert_fields(D,{'lon','lat'},'sea_ice_grid_coordinates.mat');
OUT = struct('lon',D.lon,'lat',D.lat);
save_struct(dst,OUT,OVERWRITE);
record(fid,src.seaIceCoords,dst);

%% ========================================================================
% 10. ATMOSPHERIC CONTROL DATA (Figure 10 and Figure S2)
%% ========================================================================

% ---- 10-m wind speed ----
dst = fullfile(DATA_OUT,'wind_speed_JAS_1979_2024.mat');
D = load_required(src.wind,'10-m wind speed');
if isfield(D,'V10')
    V10 = D.V10;
elseif isfield(D,'U10')
    % In the archived working file U10 contains the derived 10-m wind-speed
    % series used by the analysis. Standardize the public variable name.
    V10 = D.U10;
else
    error('Neither V10 nor U10 was found in %s',src.wind);
end
OUT = struct('V10',V10);
save_struct(dst,OUT,OVERWRITE);
record(fid,src.wind,dst);

% ---- 2-m surface air temperature ----
dst = fullfile(DATA_OUT,'surface_air_temperature_JAS_1979_2024.mat');
D = load_required(src.sat,'Surface air temperature');
if ~isfield(D,'SAT')
    error('Variable "SAT" not found in %s',src.sat);
end
OUT = struct('SAT',D.SAT);
save_struct(dst,OUT,OVERWRITE);
record(fid,src.sat,dst);

% ---- Monthly Arctic Oscillation index ----
dst = fullfile(DATA_OUT,'arctic_oscillation_monthly_1979_2024.mat');
D = load_required(src.ao,'Arctic Oscillation index');
if ~isfield(D,'AO_index')
    error('Variable "AO_index" not found in %s',src.ao);
end
OUT = struct('AO_index',D.AO_index);
save_struct(dst,OUT,OVERWRITE);
record(fid,src.ao,dst);

%% ========================================================================
% 11. SIX MONTHLY RIVER-DISCHARGE SERIES
%% ========================================================================

riverKeys = {'mackenzie','yukon','kolyma','lena','yenisei','ob'};
riverAliases = { ...
    {'mackenzie','maken'}, ...
    {'yukon'}, ...
    {'kolyma'}, ...
    {'lena'}, ...
    {'yenisei','enisey'}, ...
    {'ob'}};

searchRoots = {RIVER_DIR};

for ir = 1:numel(riverKeys)
    key = riverKeys{ir};
    manualPath = MANUAL_DISCHARGE.(key);

    if ~isempty(manualPath)
        sourceFile = manualPath;
    else
        sourceFile = find_discharge_file(searchRoots,riverAliases{ir});
    end

    if isempty(sourceFile)
        warning(['Discharge source for %s was not found automatically. ' ...
                 'Set MANUAL_DISCHARGE.%s at the top of this script.'],key,key);
        fprintf(fid,'MISSING discharge source: %s\n',key);
        continue
    end

    D = load(sourceFile);
    discharge_m3s = extract_discharge_vector(D,sourceFile);

    dst = fullfile(DATA_OUT,sprintf('discharge_%s_1972_2024.mat',key));
    OUT = struct();
    OUT.discharge_m3s = double(discharge_m3s(:));
    OUT.start_year = 1972;
    OUT.end_year   = 2024;
    OUT.temporal_resolution = 'monthly';
    save_struct(dst,OUT,OVERWRITE);
    record(fid,sourceFile,dst);
end

%% ========================================================================
% 12. SIX DRAINAGE-BASIN BOUNDARIES
%% ========================================================================

for ir = 1:numel(riverKeys)
    key = riverKeys{ir};
    sourceFile = src.basins.(key);

    if ~exist(sourceFile,'file')
        warning('Basin-boundary source not found: %s',sourceFile);
        fprintf(fid,'MISSING basin source: %s\n',sourceFile);
        continue
    end

    D = load(sourceFile);
    basin_boundary = extract_2xN_boundary(D,sourceFile);

    dst = fullfile(GEO_OUT,sprintf('basin_%s.mat',key));
    OUT = struct('basin_boundary',double(basin_boundary));
    save_struct(dst,OUT,OVERWRITE);
    record(fid,sourceFile,dst);
end

%% ========================================================================
% 13. ARCTIC MARGINAL-SEA SHAPEFILE
%% ========================================================================

if exist(src.marginalSeaShp,'file')
    copy_shapefile_set(src.marginalSeaShp,GEO_OUT,'arctic_marginal_seas',OVERWRITE);
    record(fid,src.marginalSeaShp,fullfile(GEO_OUT,'arctic_marginal_seas.shp'));
else
    warning('Marginal-sea shapefile not found: %s',src.marginalSeaShp);
    fprintf(fid,'MISSING shapefile: %s\n',src.marginalSeaShp);
end

%% ========================================================================
% 14. MAJOR-RIVER SHAPEFILE USED BY FIGURES 1 AND 7
%% ========================================================================

if ~isempty(MANUAL_MAJOR_RIVER_SHP)
else
        RIVER_DIR});
end

else
    warning(['Major-river shapefile was not found automatically. ' ...
             'Set MANUAL_MAJOR_RIVER_SHP at the top of this script.']);
    fprintf(fid,'MISSING major-river shapefile\n');

%% ========================================================================
% 15. WRITE A MACHINE-READABLE GENERATED-FILE SUMMARY
%% ========================================================================

summaryFile = fullfile(OUTPUT_ROOT,'GENERATED_FILES.txt');
write_generated_summary(summaryFile,DATA_OUT);

fprintf('\n============================================================\n');
fprintf('Finished. Generated publication inputs are in:\n%s\n',DATA_OUT);
fprintf('Preparation log:\n%s\n',logFile);
fprintf('Generated-file summary:\n%s\n',summaryFile);
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

function p = first_existing(candidates)
    p = '';
    for i = 1:numel(candidates)
        if exist(candidates{i},'file')
            p = candidates{i};
            return
        end
    end
    % Return the first candidate so a later error message remains explicit.
    if ~isempty(candidates)
        p = candidates{1};
    end
end

function ensure_dir(p)
    if ~exist(p,'dir')
        mkdir(p);
    end
end

function D = load_required(filePath,label)
    if isempty(filePath) || ~exist(filePath,'file')
        error('%s source file was not found: %s',label,filePath);
    end
    fprintf('Reading: %s\n',filePath);
    D = load(filePath);
end

function OUT = copy_fields_if_present(IN,OUT,fields)
    for i = 1:numel(fields)
        f = fields{i};
        if isfield(IN,f)
            OUT.(f) = IN.(f);
        end
    end
end

function assert_fields(S,fields,label)
    missing = {};
    for i = 1:numel(fields)
        if ~isfield(S,fields{i})
            missing{end+1} = fields{i}; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error('%s is missing required variables: %s',label,strjoin(missing,', '));
    end
end

function save_struct(filePath,S,overwrite)
    if exist(filePath,'file') && ~overwrite
        fprintf('Skipping existing file: %s\n',filePath);
        return
    end
    fprintf('Writing: %s\n',filePath);
    save(filePath,'-struct','S','-v7.3');
end

function record(fid,src,dst)
    fprintf(fid,'%s\n  -> %s\n\n',src,dst);
end

function f = find_discharge_file(searchRoots,aliases)
    f = '';
    candidates = struct('path',{},'score',{});

    % First prefer files explicitly named discharge8.mat.
    for ir = 1:numel(searchRoots)
        root = searchRoots{ir};
        if ~exist(root,'dir'), continue, end
        L = dir(fullfile(root,'**','discharge8.mat'));
        for k = 1:numel(L)
            fullp = fullfile(L(k).folder,L(k).name);
            lp = lower(fullp);
            if any(cellfun(@(a) contains(lp,lower(a)),aliases))
                candidates(end+1).path = fullp; %#ok<AGROW>
                candidates(end).score = 100;
            end
        end
    end

    % Then search more general runoff/discharge MAT files.
    if isempty(candidates)
        for ir = 1:numel(searchRoots)
            root = searchRoots{ir};
            if ~exist(root,'dir'), continue, end
            L = dir(fullfile(root,'**','*.mat'));
            for k = 1:numel(L)
                fullp = fullfile(L(k).folder,L(k).name);
                lp = lower(fullp);
                hasAlias = any(cellfun(@(a) contains(lp,lower(a)),aliases));
                hasHydro = contains(lp,'discharge') || contains(lp,'runoff');
                if hasAlias && hasHydro
                    candidates(end+1).path = fullp; %#ok<AGROW>
                    candidates(end).score = 50;
                end
            end
        end
    end

    if isempty(candidates)
        return
    end

    % Prefer the highest score, then shortest path.
    scores = [candidates.score];
    bestScore = max(scores);
    idx = find(scores == bestScore);
    if numel(idx) > 1
        lengths = arrayfun(@(q) strlength(string(candidates(q).path)),idx);
        [~,j] = min(lengths);
        idx = idx(j);
    else
        idx = idx(1);
    end
    f = candidates(idx).path;
end

function q = extract_discharge_vector(D,filePath)
    preferred = {'discharge_m3s','discharge8','discharge','runoff','Q'};
    for i = 1:numel(preferred)
        if isfield(D,preferred{i})
            v = D.(preferred{i});
            if isnumeric(v) && isvector(v)
                q = v;
                return
            end
        end
    end

    names = fieldnames(D);
    for i = 1:numel(names)
        v = D.(names{i});
        if isnumeric(v) && isvector(v) && numel(v) >= 500
            q = v;
            return
        end
    end
    error('No monthly discharge vector found in %s',filePath);
end

function B = extract_2xN_boundary(D,filePath)
    names = fieldnames(D);
    for i = 1:numel(names)
        v = D.(names{i});
        if isnumeric(v) && ndims(v)==2 && size(v,1)==2 && size(v,2)>2
            B = v;
            return
        end
    end
    error('No 2 x N basin-boundary array found in %s',filePath);
end

function copy_shapefile_set(srcShp,dstDir,newBase,overwrite)
    [srcDir,srcBase,~] = fileparts(srcShp);
    extensions = {'.shp','.shx','.dbf','.prj','.cpg','.qpj','.sbn','.sbx'};
    for i = 1:numel(extensions)
        src = fullfile(srcDir,[srcBase extensions{i}]);
        if exist(src,'file')
            dst = fullfile(dstDir,[newBase extensions{i}]);
            if exist(dst,'file') && ~overwrite
                continue
            end
            copyfile(src,dst,'f');
            fprintf('Writing: %s\n',dst);
        end
    end
end

function shp = find_major_river_shapefile(searchRoots)
    shp = '';
    candidates = {};
    for ir = 1:numel(searchRoots)
        root = searchRoots{ir};
        if ~exist(root,'dir'), continue, end
        L = dir(fullfile(root,'**','*.shp'));
        for k = 1:numel(L)
            fullp = fullfile(L(k).folder,L(k).name);
            lp = lower(fullp);
            % Prefer Natural Earth / world river products and exclude the
            % marginal-sea shapefile.
            if contains(lp,'river') && ~contains(lp,'marginal') && ~contains(lp,'basin')
                candidates{end+1} = fullp; %#ok<AGROW>
            end
        end
    end
    if isempty(candidates)
        return
    end

    % Prefer candidate names containing natural / ne_ / world.
    score = zeros(size(candidates));
    for i = 1:numel(candidates)
        lp = lower(candidates{i});
        score(i) = 1 + 4*contains(lp,'natural') + 3*contains(lp,'ne_') + 2*contains(lp,'world');
    end
    [~,idx] = max(score);
    shp = candidates{idx};
end

function write_generated_summary(filePath,dataDir)
    fid = fopen(filePath,'w');
    if fid < 0, return, end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid,'Generated public input files\n\n');
    L = dir(fullfile(dataDir,'**','*'));
    for i = 1:numel(L)
        if L(i).isdir, continue, end
        p = fullfile(L(i).folder,L(i).name);
        rel = erase(p,[dataDir filesep]);
        fprintf(fid,'%s\n',rel);
    end
end
