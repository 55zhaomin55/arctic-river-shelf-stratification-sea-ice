function S = river_shelf_systems()
%RIVER_SHELF_SYSTEMS Canonical metadata used throughout the analysis.
%
% Sampling centers are the coordinates used in the final domain-sensitivity
% analysis. Edit this file only if the final manuscript uses different centers.

S.names = { ...
    'Mackenzie-Beaufort', ...
    'Yukon-Chukchi', ...
    'Kolyma-East Siberian', ...
    'Lena-Laptev', ...
    'Yenisei-Kara', ...
    'Ob-Kara'};

S.shortNames = {'Mackenzie','Yukon','Kolyma','Lena','Yenisei','Ob'};
S.mouthLon = [-136.5, -167.0, 165.0, 130.0, 82.0, 70.0];
S.mouthLat = [  70.5,   66.0,  71.0,  74.0, 72.0, 73.0];
S.analysisRadiusKm = 450;
S.testRadiiKm = [150 300 450 600];
S.centerOffsetsDeg = [-2 -1 0 1 2];
end
