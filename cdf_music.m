% =========================================================================
% AoA-based localization using 1D MUSIC + CSI sanitization (3 and more AP triangulation)
% - Uses the SAME triangulation + CDF workflow as the FFT script
% - Replaces AoA estimator with MUSIC (endfire, 0..180, cos steering)
% - Sanitizes CSI (SpotFi-style) BEFORE AoA estimation
% - Saves error vector to ./aoa_ablation_results/
% =========================================================================
clear; clc; close all;

%%  Load data
file_path = "./wildv2/training_data_env1.h5";
n_data_to_load = 17276;
flag_print_struct = 0;
data = load_h5_structure(file_path, n_data_to_load, flag_print_struct);

%% Configuration
idxN_list = 1:1000;          % use 100 packets, max is 17276
idxF = 100;                % subcarrier, max is 234

AP_ids = [1 2 3 4 5 6];           % choose atleast 3 APs, max is 6
% assert(numel(AP_ids) == 3, 'This script assumes exactly 3 APs');

%% Constants 
c_light = physconst("LightSpeed");
ant_sep = data.ant_sep;
fc = data.center_freq;

lambda = c_light / fc;
kd = 2*pi*ant_sep / lambda;

theta_scan_deg = 0:180;
theta_scan = deg2rad(theta_scan_deg);
cos_angles = cos(theta_scan);

err_all = zeros(numel(idxN_list),1);

%% Main loop over packets
for ii = 1:numel(idxN_list)

    idxN = idxN_list(ii);
    gt_loc = data.labels(idxN,:);     % ground truth user location

    A = [];  % line normals stacked (3x2)
    b = [];  % offsets stacked (3x1)

    %% Loop over selected APs
    for jj = 1:numel(AP_ids)

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

        x = Xsan(idxF,:).';                          % [Ant x 1]
        x = x / (norm(x) + eps);

        M = length(x);
        elements = (0:M-1).';

        %% 1D MUSIC AoA (endfire, 0..180, cos steering)
        steer = exp(-1j * kd * (elements * cos_angles));   % [M x 181]

        % covariance (single snapshot, rank-1)
        R = x * x';

        [V,D] = eig(R);
        [~,id] = sort(real(diag(D)), 'descend');
        V = V(:,id);

        K = 1;                          % force for rank-1
        En = V(:,K+1:end);

        Z = En' * steer;
        Pmusic = 1 ./ (sum(abs(Z).^2,1) + eps);

        [~,imax] = max(Pmusic);
        theta_hat_deg = theta_scan_deg(imax);        % MUSIC local AoA in [0,180]

        % Convert local -> global 
        theta_local = deg2rad(180 - theta_hat_deg);
        beta_music = wrapToPi_local(psi + theta_local);

        %% Bearing line for triangulation
        % line normal for direction beta: n = [sin(beta); -cos(beta)]
        n = [sin(beta_music); -cos(beta_music)];
        A = [A; n.'];
        b = [b; n.' * c_ap(:)];
    end

    %% Least-squares triangulation
    p_est = (A.'*A) \ (A.'*b);

    %%  Localization error 
    err_all(ii) = norm(p_est.' - gt_loc);
end

%% CDF plot
figure;
cdfplot(err_all);
grid on;
set(gca, 'XScale', 'log');

xlabel('Localization error (m)');
ylabel('CDF');
title('CDF of AoA localization error (MUSIC + CSI sanitization, 3 APs)');

fprintf('\n==== MUSIC Localization Error Summary (100 packets) ====\n');
fprintf('Mean error   : %.2f m\n', mean(err_all));
fprintf('Median error : %.2f m\n', median(err_all));
fprintf('90th pct     : %.2f m\n', prctile(err_all,90));

%% Save results for for comparison later
result_dir = fullfile(pwd, 'aoa_ablation_results');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

method_name = 'MUSIC';
num_AP_used = numel(AP_ids);
fname = sprintf('loc_error_%s_%dAPs_idxN_%dto%d.mat', ...
                method_name, num_AP_used, idxN_list(1), idxN_list(end));

save(fullfile(result_dir, fname), 'err_all', 'AP_ids', 'idxN_list', 'method_name');

fprintf('Saved localization errors to:\n%s\n', fullfile(result_dir, fname));

%% Local functions
function Xsan = sanitize_csi_spotfi(X)
% SpotFi-style CSI sanitization
% X: [F x Ant] complex CSI

    % 1) Remove random packet phase using antenna 1 as reference
    ref = X(:,1);
    X1 = X .* exp(-1j*angle(ref));

    % 2) Remove linear phase trend across subcarriers (timing/CFO)
    phi = unwrap(angle(X1(:,1)));
    sc  = (1:length(phi)).';
    p   = polyfit(sc, phi, 1);
    phi_lin = polyval(p, sc);

    Xsan = X1 .* exp(-1j*phi_lin);
end

function ang = wrapToPi_local(ang)
% Wrap radians to (-pi, pi]
    ang = mod(ang + pi, 2*pi) - pi;
end
