%{
--------------------------------------------------------------------------
Το συγκεκριμένο πρόγραμμα προσομοιώνει την λειτουργία ενός FMCW ραντάρ και 
την μέτρηση απόστασης, ταχύτητας και γωνίας άφιξης του στόχου με χρήση 
Range-FFT, Doppler-FFT και Angle-FFT. Για το σήμα του ραντάρ έχει παραχθεί
ένα Radar Data Cube με όνομα beat(n,m,k). Tο n αναφέρεται στην
δειγματοληψία κάθε chirp του σήματος (fast time) και ο FFT σε αυτή την
διάσταση δίνει την beat frequency και αρα την απόσταση του στόχου. Το m
αναφέρεται στο chirp frame και ο FFT σε αυτή την διάσταση δίνει την
συχνότητα Doppler, αρα και την ταχύτητα του στόχου. Το k αναφέρεται στον
αριθμό των κεραιών του δέκτη και ο FFT σε αυτή την διάσταση δίνει την
χωρική συχνότητα u απο όπου μπορεί να υπολογιστεί το ΑοΑ. 

Κάποιες παρατηρήσεις: 
1) Το αρχικό παράδειγμα έχει γίνει για στόχο σε απόσταση 50m, ακτινική
ταχύτητα ίση με 5m/s και υπο γωνία 20 μοίρες. Μπορούν να χρησιμοποιηθούν
και διαφορετικά νούμερα, αλλά χρειάζεται προσοχή στα όρια των FFT και σε
πρακτικά όρια του ραντάρ. 

2) Έχει θεωρηθεί chirp frame με 64 chirps και ότι το ραντάρ έχει 64 κεραίες
στον δέκτη. Οι τιμές αυτές επιλέχθηκαν τυχαία, και οχι με βάση κάποιο
συγκεκριμένο ραντάρ. Αυτά μπορούν να αυξηθούν ώστε να αυξηθεί το resolution 
του FFT και να υπάρχει μικρότερο σφάλμα.
---------------------------------------------------------------------------
Ιωάννης Βεκίλογλου
Πανεπιστήμιο Δυτικής Αττικής
28/05/26
---------------------------------------------------------------------------
%}

clear; clc; close all;
%% ============================================================
%                    FMCW RADAR PARAMETERS
% ============================================================

c  = 3e8;
fc = 5.3e9;
lambda = c / fc;

BW = 600e6;
Tc = 1e-3;
S  = BW / Tc;

fs = 500e3;
N  = round(fs * Tc);

M  = 64;     % Αριθμός chirps
K  = 64;     % Κεραίες δέκτη

chirp_idx = 0:M-1;
ant_idx   = 0:K-1;

%% ============================================================
%                          TARGET
% ============================================================

R_target     = 50;
v_target     = 5; % Ακτινική ταχύτητα
theta_target = 20;

tau = 2 * R_target / c;
f_D = -(2 * v_target) / lambda; % Προσοχή στον τύπο της συχνότητας doppler

theta_rad = deg2rad(theta_target);
d = lambda/2;
dphi = (2*pi*d*sin(theta_rad))/lambda;

fprintf('Beat Frequency: %.2f kHz\n', (2*S*R_target/c)/1e3);
fprintf('Doppler Frequency: %.2f Hz\n', f_D);
fprintf('Spatial phase step: %.4f rad\n', dphi);

%% ============================================================
%                    FAST TIME BASE
% ============================================================

t_fast = (0:N-1)/fs;

%% ============================================================
%                    TX SIGNAL
% ============================================================

phi_tx = 2*pi*(fc*t_fast + 0.5*S*t_fast.^2);
tx = exp(1j * phi_tx); % Σήμα Chirp στον πομπό του ραντάρ.

%% ============================================================
%                    SIGNAL CUBE (N x M x K)
% ============================================================

beat = zeros(N, M, K);

for m = 1:M

    doppler = exp(1j * 2*pi * f_D * (m-1) * Tc);

    t_delayed = t_fast - tau;

    rx = exp(1j * 2*pi*(fc*t_delayed + 0.5*S*t_delayed.^2));

    for k = 1:K

        spatial_phase = exp(1j * (k-1) * dphi);

        beat(:,m,k) = tx .* conj(rx) * doppler * spatial_phase; 

    end
end


%% ============================================================
%                        RANGE FFT
% ============================================================

range_fft = fft(beat, N, 1);
range_fft = range_fft(1:N/2,:,:);

range_mag = abs(range_fft);
range_mag = range_mag / max(range_mag(:));

f_range = (0:N/2-1)*(fs/N); % Στον FFT ισχύει: f[k] = k*(fs/N)
range_axis = (f_range * c)/(2*S);

%% ============================================================
%                        DOPPLER FFT
% ============================================================

RD = fftshift(fft(range_fft, M, 2), 2);

RD_mag = abs(RD);
RD_mag = RD_mag / max(RD_mag(:));

RD_dB = 20*log10(RD_mag + 1e-12);

PRF = 1/Tc;
fd_axis = (-M/2:M/2-1)*(PRF/M);
vel_axis = -(lambda * fd_axis)/2;

%% ============================================================
%                        ANGLE FFT
% ============================================================

RA = fftshift(fft(RD, K, 3), 3);

RA_mag = abs(RA);
RA_mag = RA_mag / max(RA_mag(:));
RA_dB = 20*log10(RA_mag + 1e-12);

u_axis = (-K/2:K/2-1)/K;

sin_theta = (u_axis * lambda)/d;
sin_theta = max(min(sin_theta,1),-1);
theta_axis = asind(sin_theta);

%% ============================================================
%                  PEAK DETECTION (3D)
% ============================================================

[~, idx] = max(RA_mag(:));
[row, col, ant] = ind2sub(size(RA_mag), idx);

R_est     = range_axis(row);
v_est     = vel_axis(col);
theta_est = theta_axis(ant);


%% ============================================================
%              CLEAN FMCW CHIRP SPECTROGRAM
% ============================================================

fs_vis = 2e9;               % High fs just for display (must be > 2*BW = 1.2 GHz)
Tc_vis = Tc;                  % Same chirp duration (1 ms)
N_vis  = round(fs_vis * Tc_vis);   % Number of samples

t_vis  = (0:N_vis-1) / fs_vis;    % Fast-time axis at high rate


phi_vis = 2*pi * (0.5 * S * t_vis.^2);   % S = BW/Tc 
tx_vis  = exp(1j * phi_vis);


%{
Εδω κάνω plot μόνο το chirp (δλδ απο 0 - 600MHz), 
για να φανεί όλο το σήμα κάνε uncomment out τις γραμμές απο κάτω
%} 

figure; 
spectrogram(tx_vis, ...
    hamming(512), ...          
    480, ...                   
    2048, ...                  
    fs_vis, ...                
'yaxis');
title('FMCW Chirp Spectrogram (BW = 600 MHz, T_c = 1 ms)');
colormap turbo;

% [s, f, t] = spectrogram(tx_vis, ...   % tx_vis uses only 0.5*S*t.^2 (no fc)
%      hamming(512), 480, 2048, fs_vis);
% 
%  figure;
%  imagesc(t*1e3, (f + fc)/1e9, 20*log10(abs(s) + 1e-12));
%  axis xy;
%  colormap turbo;
%  colorbar;
% 
% xlabel('Time (ms)');
% ylabel('Frequency (GHz)');
% title('FMCW Chirp Spectrogram (f_c = 5.3 GHz, BW = 600 MHz)');
%% ============================================================
%           RANGE FFT (FREQUENCY DOMAIN)
% ============================================================

figure;
plot(f_range/1e3, range_mag(:,1,1), 'LineWidth',1.5);
grid on;
xlabel('Beat Frequency (kHz)');
ylabel('Magnitude');
title('Range FFT (Frequency Domain)');
hold on;

[~, peak_r] = max(range_mag(:,1,1));
plot(f_range(peak_r)/1e3, range_mag(peak_r,1,1), ...
    'ro','MarkerFaceColor','r');

text(f_range(peak_r)/1e3, range_mag(peak_r,1,1), ...
    sprintf(' %.1f kHz', f_range(peak_r)/1e3));

%% ============================================================
%              RANGE FFT (RANGE DOMAIN)
% ============================================================

figure;
plot(range_axis, range_mag(:,1,1), 'LineWidth',1.5);
grid on;
xlabel('Range (m)');
ylabel('Magnitude');
title('Range FFT (Range Domain)');
hold on;

plot(R_est, range_mag(peak_r,1,1), ...
    'ro','MarkerFaceColor','r');

text(R_est, range_mag(peak_r,1,1), ...
    sprintf(' %.2f m', R_est));

%% ============================================================
%           DOPPLER FFT (FREQUENCY DOMAIN)
% ============================================================

doppler_slice = squeeze(RD_mag(peak_r,:,1));
[~, peak_d] = max(doppler_slice);

figure;
plot(fd_axis, doppler_slice, 'LineWidth',1.5);
grid on;
xlabel('Doppler Frequency (Hz)');
ylabel('Magnitude');
title('Doppler FFT');
hold on;

plot(fd_axis(peak_d), doppler_slice(peak_d), ...
    'ro','MarkerFaceColor','r');

text(fd_axis(peak_d), doppler_slice(peak_d), ...
    sprintf(' %.1f Hz', fd_axis(peak_d)));

%% ============================================================
%           DOPPLER FFT (VELOCITY DOMAIN)
% ============================================================

figure;
plot(vel_axis, doppler_slice, 'LineWidth',1.5);
grid on;

xlabel('Velocity (m/s)');
ylabel('Magnitude');
title('Doppler FFT (Velocity Domain)');

hold on;
plot(v_est, doppler_slice(peak_d), ...
    'ro','MarkerFaceColor','r');

text(v_est, doppler_slice(peak_d), ...
    sprintf(' %.2f m/s', v_est));

%% ============================================================
%           ANGLE FFT 
% ============================================================

angle_slice = squeeze(RA_mag(peak_r, peak_d, :));

figure;
plot(theta_axis, angle_slice, 'LineWidth',1.8);
grid on;
xlabel('Angle (deg)');
ylabel('Magnitude');
title('Angle FFT');
hold on;

plot(theta_est, angle_slice(ant), ...
    'ro','MarkerFaceColor','r');

text(theta_est, angle_slice(ant), ...
    sprintf(' %.2f°', theta_est));

%% ============================================================
%                 RANGE–ANGLE MAP
% ============================================================

range_angle = squeeze(max(RA_dB, [], 2)); 

figure;
imagesc(theta_axis, range_axis, range_angle);
axis xy;
colormap turbo;
colorbar;
caxis([-40 0]);

xlabel('Angle (deg)');
ylabel('Range (m)');
title('Range–Angle Map');

% Αν θέλεις marker κάνε uncomment out τις παρακάτω γραμμές:

% hold on;
% plot(theta_est, R_est, 'ro','MarkerFaceColor','r');
% 
% text(theta_est, R_est, ...
%     sprintf(' (%.1f m, %.1f°)', R_est, theta_est), ...
%     'Color','w','FontWeight','bold');

%% ============================================================
%                 RANGE–DOPPLER MAP
% ============================================================

range_doppler = squeeze(max(RD_dB, [], 3)); 

figure;

imagesc(vel_axis, range_axis, range_doppler);

axis xy;

colormap turbo;
colorbar;

caxis([-40 0]);

xlabel('Velocity (m/s)');
ylabel('Range (m)');

title('Range–Doppler Map');

% Αν θέλεις marker κάνε uncomment out τις παρακάτω γραμμές:

% hold on;
% 
% plot(v_est, R_est, ...
%     'ro','MarkerFaceColor','r');
% 
% text(v_est, R_est, ...
%     sprintf(' (%.1f m, %.1f m/s)', ...
%     R_est, v_est), ...
%     'Color','w','FontWeight','bold');

%% ============================================================
%                 DOPPLER–ANGLE MAP
% ============================================================

doppler_angle = squeeze(max(RA_dB, [], 1)); 

figure;

imagesc(theta_axis, vel_axis, doppler_angle);

axis xy;

colormap turbo;
colorbar;

caxis([-40 0]);

xlabel('Angle (deg)');
ylabel('Velocity (m/s)');

title('Doppler–Angle Map');

% Αν θέλεις marker κάνε uncomment out τις παρακάτω γραμμές:

% hold on;
% 
% plot(theta_est, v_est, ...
%     'ro','MarkerFaceColor','r');
% 
% text(theta_est, v_est, ...
%     sprintf(' (%.1f°, %.1f m/s)', ...
%     theta_est, v_est), ...
%     'Color','w','FontWeight','bold');

%% ============================================================
%                    FINAL RESULTS
% ============================================================
range_error  = abs(R_target - R_est);
vel_error    = abs(v_target - v_est);
angle_error  = abs(theta_target - theta_est);

fprintf('\n====================================\n');
fprintf('         ESTIMATION RESULTS\n');
fprintf('====================================\n');

fprintf('RANGE:\n');
fprintf('  True      : %.2f m\n', R_target);
fprintf('  Estimated : %.2f m\n', R_est);
fprintf('  Error     : %.4f m\n\n', range_error);

fprintf('VELOCITY:\n');
fprintf('  True      : %.2f m/s\n', v_target);
fprintf('  Estimated : %.2f m/s\n', v_est);
fprintf('  Error     : %.4f m/s\n\n', vel_error);

fprintf('ANGLE:\n');
fprintf('  True      : %.2f deg\n', theta_target);
fprintf('  Estimated : %.2f deg\n', theta_est);
fprintf('  Error     : %.4f deg\n', angle_error);

fprintf('====================================\n');
