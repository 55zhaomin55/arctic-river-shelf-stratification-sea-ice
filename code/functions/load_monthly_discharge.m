function discharge = load_monthly_discharge(dataDir)
%LOAD_MONTHLY_DISCHARGE Load six monthly ArcticGRO discharge series.
%
% Returns a 6 x N matrix in the system order used by the manuscript.

files = { ...
    'discharge_mackenzie_1972_2024.mat', ...
    'discharge_yukon_1972_2024.mat', ...
    'discharge_kolyma_1972_2024.mat', ...
    'discharge_lena_1972_2024.mat', ...
    'discharge_yenisei_1972_2024.mat', ...
    'discharge_ob_1972_2024.mat'};

discharge = [];
for i = 1:numel(files)
    f = fullfile(dataDir,files{i});
    if ~exist(f,'file')
        error('Missing discharge file: %s',f);
    end
    D = load(f);
    if isfield(D,'discharge8')
        q = D.discharge8;
    elseif isfield(D,'discharge_m3s')
        q = D.discharge_m3s;
    else
        names = fieldnames(D);
        q = [];
        for k = 1:numel(names)
            v = D.(names{k});
            if isnumeric(v) && isvector(v)
                q = v;
                break
            end
        end
        if isempty(q)
            error('No numeric discharge vector found in %s',f);
        end
    end
    discharge(i,1:numel(q)) = double(q(:)); %#ok<AGROW>
end
end
