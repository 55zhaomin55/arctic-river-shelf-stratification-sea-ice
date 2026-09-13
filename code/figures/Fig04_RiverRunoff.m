%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% =========================
% 1. Load monthly discharge
% =========================
dischargeAll = load_monthly_discharge(cfg.dataDir); % 6 x 636, Jan 1972-Dec 2024
discharge1 = dischargeAll(:,85:end)';                    % 552 x 6, Jan 1979-Dec 2024

%% =========================
% 2. Basic information
% =========================
years = 1979:2024;
nYear = length(years);

riverNames = {'Mackenzie','Yukon','Kolyma','Lena','Yenisei','Ob'};

% Days per month (May-September is unaffected by leap-year February)
days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];

% Target season: May-September
target_months = 5:9;   % May–Sep
target_month_labels = 'May–Sep';

%% =========================
% 3. Calculate cumulative May-September runoff volume (km^3)
%    A year is missing if any target month is missing.
% =========================
Q_MaySep_km3 = nan(nYear, 6);

for y = 1:nYear
    for b = 1:6

        idx_months = (y-1)*12 + target_months;
        q_months = discharge1(idx_months, b);

        if any(isnan(q_months))
            Q_MaySep_km3(y, b) = NaN;
            continue
        end

        total_volume = 0;
        for k = 1:length(target_months)
            m = target_months(k);
            sec_this_month = days_in_month(m) * 24 * 3600;
            total_volume = total_volume + q_months(k) * sec_this_month;
        end

        Q_MaySep_km3(y, b) = total_volume / 1e9;
    end
end

%% =========================
% 4. Report missing years
% =========================
disp('===== Missing years for May–Sep cumulative runoff =====')
for i = 1:6
    fprintf('%s missing years:\n', riverNames{i});
    disp(years(isnan(Q_MaySep_km3(:, i)))');
end

%% =========================
% 5. Plot
% =========================
figure(1); clf
set(gcf, 'Units', 'centimeters', 'Position', [2 2 17 17], 'Color', 'w');

% ====== subplot positions ======
pos1 = [0.05, 0.67, 0.40, 0.26];
pos2 = [0.52, 0.67, 0.40, 0.26];
pos3 = [0.05, 0.36, 0.40, 0.26];
pos4 = [0.52, 0.36, 0.40, 0.26];
pos5 = [0.05, 0.05, 0.40, 0.26];
pos6 = [0.52, 0.05, 0.40, 0.26];
pos = [pos1; pos2; pos3; pos4; pos5; pos6];

linewidth_raw    = 0.8;
linewidth_smooth = 1.6;
font = 6;

for i = 1:6

    subplot('Position', pos(i,:)); hold on;

    y = Q_MaySep_km3(:, i);

    % 5-year moving mean
    y_smooth = movmean(y, 5, 'omitnan', 'Endpoints', 'shrink');

    h1 = plot(years, y, '-', ...
        'Color', [0.00 0.55 0.80], ...
        'LineWidth', linewidth_raw);

    h2 = plot(years, y_smooth, '-', ...
        'Color', [0.00 0.30 0.90], ...
        'LineWidth', linewidth_smooth);

    idx = isfinite(y);
    if sum(idx) >= 3
        p = polyfit(years(idx), y(idx), 1);
        y_fit = polyval(p, years);
        h3 = plot(years, y_fit, '--', ...
            'Color', [0.3 0.3 0.3], ...
            'LineWidth', 0.9);

        [~, P] = corrcoef(years(idx), y(idx));
        pval = P(1,2);

        slope = p(1);   % km^3 / year

        if pval < 0.001
            stars = '***';
        elseif pval < 0.01
            stars = '**';
        elseif pval < 0.05
            stars = '*';
        else
            stars = '';
        end

        trend_str = sprintf('Trend = %.2f km^3 yr^{-1} %s', slope, stars);
    else
        h3 = plot(nan, nan, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.9);
        trend_str = 'Trend = NaN';
    end

    xlim([1975 2027]);
    xticks([1980 1990 2000 2010 2020]);

    y_valid = y(isfinite(y));
    if isempty(y_valid)
        ymin = 0;
        ymax = 100;
    else
        ymin = floor(min(y_valid)*0.9/10)*10;
        ymax = ceil(max(y_valid)*1.1/10)*10;
        if ymin == ymax
            ymin = ymin - 10;
            ymax = ymax + 10;
        end
    end
    ylim([ymin ymax]);

    text(1976, ymax - 0.10*(ymax-ymin), riverNames{i}, ...
        'FontSize', font+5, ...
        'FontName', 'Times New Roman', ...
        'FontWeight', 'bold');

    text(1996, ymin + 0.08*(ymax-ymin), trend_str, ...
        'FontSize', font+1, ...
        'FontName', 'Times New Roman', ...
        'Color', [0.2 0.2 0.2], ...
        'FontWeight', 'bold');

    if ismember(i, [1 3 5])
        ylabel('Cumulative runoff (km^3)', ...
            'FontName', 'Times New Roman', ...
            'FontSize', font+1, ...
            'FontWeight', 'bold');
    end

    if ismember(i, [5 6])
        xlabel('Year', ...
            'FontName', 'Times New Roman', ...
            'FontSize', font+1, ...
            'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end

    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', font, ...
        'FontWeight', 'bold', ...
        'LineWidth', 0.8, ...
        'Box', 'on', ...
        'TickDir', 'out', ...
        'XColor', 'k', ...
        'YColor', 'k');

    grid on; box on;

    if i == 1
        h = legend([h1, h2, h3], ...
            {['Annual ', target_month_labels, ' cumulative runoff'], ...
             '5-year moving mean', ...
             'Linear trend'}, ...
            'Location', 'east', ...
            'Box', 'off', ...
            'FontSize', 7);
        h.ItemTokenSize = [5, 5];
    end
end

%% =========================
% =========================
set(gcf,'InvertHardcopy','off');
exportgraphics(gcf,fullfile(cfg.outputDir,'Fig04_RiverRunoff.png'),'Resolution',600);
