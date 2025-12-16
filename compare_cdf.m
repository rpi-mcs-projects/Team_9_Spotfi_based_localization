clear; clc; close all;

%% Result directory 
result_dir = fullfile(pwd, 'aoa_ablation_results');

%% Files to compare 
files = {
    % 'loc_error_FFT_3APs_idxN_1to55.mat'
    % 'loc_error_FFT_4APs_idxN_1to55.mat'
    % 'loc_error_FFT_5APs_idxN_1to55.mat'
    'loc_error_FFT_6APs_idxN_1to55.mat'
    'loc_error_FFT_6APs_idxN_1to100.mat'
    'loc_error_FFT_6APs_idxN_1to500.mat'
    % 'loc_error_MUSIC_3APs_idxN_1to55.mat'
    % 'loc_error_MUSIC_4APs_idxN_1to55.mat'
    % 'loc_error_MUSIC_5APs_idxN_1to55.mat'
    % 'loc_error_MUSIC_6APs_idxN_1to55.mat'
    
};

%% Plot CDFs
figure; hold on; grid on;

legend_entries = cell(numel(files),1);

for k = 1:numel(files)

    data = load(fullfile(result_dir, files{k}));

    % Plot CDF
    h = cdfplot(data.err_all);
    h.LineWidth = 2;

    % Create clean legend label from filename
    % Example: loc_error_FFT_3APs_idxN_1to55.mat → FFT (3 APs)
    tokens = regexp(files{k}, 'loc_error_(.*)_([0-9]+)APs', 'tokens');
    if ~isempty(tokens)
        method = tokens{1}{1};
        nAPs   = tokens{1}{2};
        legend_entries{k} = sprintf('%s (%s APs)', method, nAPs);
    else
        legend_entries{k} = files{k};
    end
end

%% Improve plot aesthetics
set(gcf, 'Color', 'w');                 % white background
set(gcf, 'Position', [100 100 700 500]);% larger figure size

ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;

set(gca, 'XScale', 'log');

xlabel('Localization error (m)', ...
       'FontSize', 16, 'FontWeight', 'bold');
ylabel('CDF', ...
       'FontSize', 16, 'FontWeight', 'bold');

title('CDF comparison: FFT vs MUSIC (3 APs)', ...
      'FontSize', 16, 'FontWeight', 'bold');

legend(legend_entries, ...
       'Location', 'southeast', ...
       'FontSize', 13, ...
       'Box', 'off');

grid on;
