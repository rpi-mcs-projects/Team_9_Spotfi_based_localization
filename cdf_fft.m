% =========================================================================
% AoA-based localization using FFT + CSI sanitization
% - Sanitizes CSI (SpotFi-style) BEFORE AoA estimation
% - Saves error vector to ./aoa_ablation_results/
% =========================================================================
clear; clc; close all;

%% Load data
file_path = "./wildv2/training_data_env1.h5";
n_data_to_load = 17276;
flag_print_struct = 0;
data = load_h5_structure(file_path, n_data_to_load, flag_print_struct);

%% Result directory 
result_dir = fullfile(pwd, 'aoa_ablation_results');

if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

%% Configuration 
idxN_list = 1:1000;          % use 100 packets, max is 17276
idxF = 100;                % subcarrier, max is 234

AP_ids = [1 2 3 4 5 6];           % choose atleast any 3 APs you want (max:6)
numAP = numel(AP_ids);


%%  Constants
c_light = physconst("LightSpeed");
ant_sep = data.ant_sep;
fc = data.center_freq;

lambda = c_light / fc;
kd = 2*pi*ant_sep / lambda;

Nfft = 1024;               % FFT size
err_all = zeros(numel(idxN_list),1);

%% Main loop over packets
for ii = 1:numel(idxN_list)

    idxN = idxN_list(ii);
    gt_loc = data.labels(idxN,:);     % ground truth user location

    A = [];    % normal vectors for lines
    b = [];    % offsets

    %% Loop over selected APs
    for jj = 1:numAP

        idxAP = AP_ids(jj);

        %% AP geometry
        ap_loc = squeeze(data.ap_locs(idxAP,:,:));   % [4x2]
        c_ap = mean(ap_loc,1);                       % centroid

        % array axis (ant1 -> ant4)
        v = ap_loc(end,:) - ap_loc(1,:);
        psi = atan2(v(2), v(1));                     % global axis angle

        %% CSI sanitization
        Xraw = squeeze(data.cmplx_channel(idxN,:, :, idxAP));  % [F x Ant]
        Xsan = sanitize_csi_spotfi(Xraw);

        x = Xsan(idxF,:).';
        x = x / norm(x);

        %% FFT AoA (endfire, 0..180)
        Xf = fftshift(fft(x, Nfft));
        u = linspace(-1,1,Nfft);          % u = cos(theta)
        theta_fft = acosd(u);             % [0..180]

        Pfft = abs(Xf).^2;
        [~,imax] = max(Pfft);
        theta_hat = theta_fft(imax);      % local AoA (deg)

        % convert local -> global (your convention)
        theta_local = deg2rad(180 - theta_hat);
        beta_fft = wrapToPi_local(psi + theta_local);

        %% Line equation for triangulation
        % line normal: n = [sin(beta); -cos(beta)]
        n = [sin(beta_fft); -cos(beta_fft)];

        A = [A; n.'];
        b = [b; n.' * c_ap(:)];
    end

    %% Least-squares triangulation
    % Solve A * p = b
    p_est = (A.'*A) \ (A.'*b);

    %% Localization error
    err_all(ii) = norm(p_est.' - gt_loc);

end

%% Save error results 
method_name = 'FFT';
num_AP_used = numel(AP_ids);

% Create filename and save it for CDF comparison later
fname = sprintf('loc_error_%s_%dAPs_idxN_%dto%d.mat', ...
                method_name, num_AP_used, idxN_list(1), idxN_list(end));

save(fullfile(result_dir, fname), ...
     'err_all', 'AP_ids', 'idxN_list', 'method_name');

fprintf('Saved localization errors to:\n%s\n', fullfile(result_dir, fname));


%% CDF plot
figure;
cdfplot(err_all);
grid on;

set(gca, 'XScale', 'log'); 

xlabel('Localization error (m)');
ylabel('CDF');
title('CDF of AoA-based localization error (FFT + CSI sanitization)');

fprintf('\n==== Localization Error Summary (100 packets) ====\n');
fprintf('Mean error   : %.2f m\n', mean(err_all));
fprintf('Median error : %.2f m\n', median(err_all));
fprintf('90th pct     : %.2f m\n', prctile(err_all,90));

%% Local functions 
function Xsan = sanitize_csi_spotfi(X)
% SpotFi-style CSI sanitization
% X: [F x Ant]

    % remove random packet phase
    ref = X(:,1);
    X1 = X .* exp(-1j*angle(ref));

    % remove linear phase slope across subcarriers
    phi = unwrap(angle(X1(:,1)));
    sc  = (1:length(phi)).';
    p   = polyfit(sc, phi, 1);
    phi_lin = polyval(p, sc);

    Xsan = X1 .* exp(-1j*phi_lin);
end

function ang = wrapToPi_local(ang)
% wrap radians to (-pi, pi]
    ang = mod(ang + pi, 2*pi) - pi;
end
