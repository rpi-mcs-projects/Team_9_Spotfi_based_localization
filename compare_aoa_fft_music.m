% =========================================================================
% Compare AoA estimation of FFT Vs MUSIC (both single dimension AoA)
% Ground truth angle wrt to each AP towards a specific location is shown as
% well
% prints the mean error from each method.
% =========================================================================
clear; close all; clc;

%% Load data 
file_path = "./wildv2/training_data_env1.h5";
n_data_to_load = 17276;
flag_print_struct = 1;
data = load_h5_structure(file_path, n_data_to_load, flag_print_struct);

%% Choose snapshot 
idxN = 60;          % packet index
idxF = 100;         % subcarrier index (after sanitization)

gt_loc = data.labels(idxN,:);      % user location [1x2]
numAP  = size(data.ap_locs,1);

%% Constants 
c_light = physconst("LightSpeed");
ant_sep = data.ant_sep;
fc      = data.center_freq;

lambda = c_light / fc;
kd = 2*pi*ant_sep / lambda;

theta_scan_deg = 0:180;
theta_scan = deg2rad(theta_scan_deg);
cos_angles = cos(theta_scan);

%% Plot setup 
figure; hold on; axis equal; grid on;
xlabel('x'); ylabel('y');
title('All APs: Geometry (green) | MUSIC (cyan) | FFT (magenta)');

plot(gt_loc(1), gt_loc(2), 'md', 'MarkerFaceColor','m', 'MarkerSize',8);
text(gt_loc(1), gt_loc(2), ' USER', 'Color','m', 'FontWeight','bold');

L = 10;   % ray length

% colors for geom(ground truth), music, and fft
col_geom  = [0 0.7 0];
col_music = [0 0.8 0.8];
col_fft   = [0.8 0 0.8];

beta_err_music = [];
beta_err_fft   = [];

%% Loop over APs 
for idxAP = 1:numAP

    %% Geometry
    ap_loc = squeeze(data.ap_locs(idxAP,:,:));   % [4x2]
    plot(ap_loc(:,1), ap_loc(:,2), 'ko-', 'LineWidth',1.5,'MarkerFaceColor','k');

    c = mean(ap_loc,1);
    plot(c(1),c(2),'bs','MarkerFaceColor','b');
    text(c(1),c(2),sprintf(' AP%d',idxAP),'Color','b');

    % Array axis (ant1 -> ant4)
    v = ap_loc(end,:) - ap_loc(1,:);
    psi = atan2(v(2),v(1));   % global

    % GT global bearing
    beta_gt = atan2(gt_loc(2)-c(2), gt_loc(1)-c(1));
    p_gt = c + L*[cos(beta_gt), sin(beta_gt)];
    plot([c(1) p_gt(1)], [c(2) p_gt(2)], '--', 'Color',col_geom,'LineWidth',2);

    %% CSI sanitization 
    Xraw = squeeze(data.cmplx_channel(idxN,:, :, idxAP));   % [F x Ant]
    Xsan = sanitize_csi_spotfi(Xraw);
    x = Xsan(idxF,:).';
    x = x / norm(x);

    M = length(x);
    elements = (0:M-1).';

    %% 1D MUSIC 
    steer = exp(-1j * kd * (elements * cos_angles));
    S = x*x';

    [V,D] = eig(S);
    [~,id] = sort(real(diag(D)),'descend');
    V = V(:,id);

    K = 1;
    En = V(:,K+1:end);

    Z = En' * steer;
    Pmusic = 1 ./ (sum(abs(Z).^2,1) + eps);

    [~,imax] = max(Pmusic);
    theta_music_hat = theta_scan_deg(imax);

    theta_music = deg2rad(180 - theta_music_hat);
    beta_music  = wrapToPi_local(psi + theta_music);

    p_mu = c + L*[cos(beta_music), sin(beta_music)];
    plot([c(1) p_mu(1)], [c(2) p_mu(2)], '--', 'Color',col_music,'LineWidth',2);

    beta_err_music(end+1) = abs(rad2deg(wrapToPi_local(beta_music - beta_gt)));

    %% 1D FFT 
    Nfft = 1024;
    Xf = fftshift(fft(x, Nfft));
    u = linspace(-1,1,Nfft);              % u = cos(theta)
    theta_fft = acosd(u);                 % 0..180
    Pfft = abs(Xf).^2;
    [~,ifk] = max(Pfft);
    theta_fft_hat = theta_fft(ifk);

    theta_fft_loc = deg2rad(180 - theta_fft_hat);
    beta_fft = wrapToPi_local(psi + theta_fft_loc);

    p_fft = c + L*[cos(beta_fft), sin(beta_fft)];
    plot([c(1) p_fft(1)], [c(2) p_fft(2)], '--', 'Color',col_fft,'LineWidth',2);

    beta_err_fft(end+1) = abs(rad2deg(wrapToPi_local(beta_fft - beta_gt)));

    fprintf('AP %d | beta_GT=%+.1f | MUSIC=%+.1f | FFT=%+.1f\n', ...
        idxAP, rad2deg(beta_gt), rad2deg(beta_music), rad2deg(beta_fft));

end

legend({'USER','Antennas','Centroid','GT AoA','MUSIC AoA','FFT AoA'}, 'Location','best');

%% Error summary
fprintf('\n=== Average Absolute Beta Error ===\n');
fprintf('MUSIC: %.2f deg\n', mean(beta_err_music));
fprintf('FFT  : %.2f deg\n', mean(beta_err_fft));

%% Local functions
function Xsan = sanitize_csi_spotfi(X)
    ref = X(:,1);
    X1 = X .* exp(-1j*angle(ref));

    phi = unwrap(angle(X1(:,1)));
    sc  = (1:length(phi)).';
    p   = polyfit(sc,phi,1);
    phi_lin = polyval(p,sc);

    Xsan = X1 .* exp(-1j*phi_lin);
end

function ang = wrapToPi_local(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end
