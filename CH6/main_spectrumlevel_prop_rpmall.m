function main_spectrumlevel_prop_rpmall
save_data = 1;
include_background = 1;
addpath('O:\MATLAB Signal Processing Files');

datafile1 = 'O:\V-Tunnel 12-12 Wing-Prop\Data MAT\Meas1.mat';
datafile2 = 'O:\V-Tunnel 12-12 Wing-Prop\Data MAT\Meas3.mat';
datafile3 = 'O:\V-Tunnel 12-12 Wing-Prop\Data MAT\Meas5.mat';

datafile4 = 'O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\Meas1.mat';
datafile5 = 'O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\Meas3.mat';
datafile6 = 'O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\Meas5.mat';

mainresultpath = 'O:\PhD Thesis\RESULTS\CH6\PROPWING';
backgroundfile = 'O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\0.Background.mat';

load('O:\V-Tunnel 12-12 Wing-Prop\mic_poses_optim.mat');
mic_config = mic_poses.'; clear mic_poses;
mic_config(:,3) = 0.02;

mic_sel = 41:45;
Fs = 50e3; % sample frequency
t_b = 1;
fband = 5;

if include_background&&exist(backgroundfile ,'file')
    fprintf('\tIncluding background data...\n');
    load(backgroundfile);
    data = data(:,mic_sel);
    
    [N, n_mic] = size(data);
    t = N/50e3; % sec
    N_b = t_b*Fs;

    n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
    fprintf('\tUsing %d blocks\n', n_blocks);
    df = Fs/N_b;

    psdxb = zeros(floor(N_b/2), n_mic);
    freqb = (0:floor(N_b/2)-1)*df;
    reverseStr = '';
    for B = 1:n_blocks
        msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        xdft = zeros(N_b, n_mic);
        for M = 1:n_mic
            % fft with hanning window to all microphones, index according to 50% overlap
            xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
                hann(N_b));
        end
        psdxb = psdxb + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

    end
    psdxb = (8/3)*psdxb/n_blocks; % Hanning AMPLITUDE CORRECTION factor
    psdxb(2:end-1,:) = 2*psdxb(2:end-1,:); % Single-side

    fprintf([reverseStr, '\tAll blocks evaluated!\n']);
    
%     mean_logpsdb = 10*log10(mean(psdxb,2));
    [slb, fb] = PSDToSpecLev(mean(psdxb,2), fband, Fs);
end

fprintf(['\tStart spectrum level calculations for signals...\n']);
clear data psdxb xdft;

load(datafile1);
data = data(:,mic_sel);

[N, n_mic] = size(data);
t = N/50e3; % sec
N_b = t_b*Fs;

n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
fprintf('\tUsing %d blocks\n', n_blocks);
df = Fs/N_b;

psdx1 = zeros(floor(N_b/2), n_mic);
freq = (0:floor(N_b/2)-1)*df;
reverseStr = [];
for B = 1:n_blocks
    msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    xdft = zeros(N_b, n_mic);
    for M = 1:n_mic
        % fft with hanning window to all microphones, index according to 50% overlap
        xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx1 = psdx1 + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx1 = (8/3)*psdx1/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx1(2:end-1,:) = 2*psdx1(2:end-1,:); % Single-side

fprintf([reverseStr, '\tAll blocks evaluated!\n']);
    
% mean_logpsd1 = 10*log10(mean(psdx1,2));
[sl1, f1] = PSDToSpecLev(mean(psdx1,2), fband, Fs);

clear data psdx1 xdft;

load(datafile2);
data = data(:,mic_sel);

[N, n_mic] = size(data);
t = N/50e3; % sec
N_b = t_b*Fs;

n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
fprintf('\tUsing %d blocks\n', n_blocks);
df = Fs/N_b;

psdx2 = zeros(floor(N_b/2), n_mic);
freq = (0:floor(N_b/2)-1)*df;
reverseStr = [];
for B = 1:n_blocks
    msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    xdft = zeros(N_b, n_mic);
    for M = 1:n_mic
        % fft with hanning window to all microphones, index according to 50% overlap
        xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx2 = psdx2 + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx2 = (8/3)*psdx2/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx2(2:end-1,:) = 2*psdx2(2:end-1,:); % Single-side

fprintf([reverseStr, '\tAll blocks evaluated!\n']);
    
% mean_logpsd2 = 10*log10(mean(psdx2,2));
[sl2, f2] = PSDToSpecLev(mean(psdx2,2), fband, Fs);

clear data psdx2 xdft;

load(datafile3);
data = data(:,mic_sel);

[N, n_mic] = size(data);
t = N/50e3; % sec
N_b = t_b*Fs;

n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
fprintf('\tUsing %d blocks\n', n_blocks);
df = Fs/N_b;

psdx3 = zeros(floor(N_b/2), n_mic);
freq = (0:floor(N_b/2)-1)*df;
reverseStr = [];
for B = 1:n_blocks
    msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    xdft = zeros(N_b, n_mic);
    for M = 1:n_mic
        % fft with hanning window to all microphones, index according to 50% overlap
        xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx3 = psdx3 + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx3 = (8/3)*psdx3/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx3(2:end-1,:) = 2*psdx3(2:end-1,:); % Single-side

fprintf([reverseStr, '\tAll blocks evaluated!\n']);
    
[sl3, f3] = PSDToSpecLev(mean(psdx3,2), fband, Fs);

clear data psdx3 xdft;

% CHECK MORE

load(datafile4);
data = data(:,mic_sel);
[N, n_mic] = size(data);
t = N/50e3; % sec
N_b = t_b*Fs;
n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
fprintf('\tUsing %d blocks\n', n_blocks);
df = Fs/N_b;
psdx = zeros(floor(N_b/2), n_mic);
freq = (0:floor(N_b/2)-1)*df;
reverseStr = [];
for B = 1:n_blocks
    msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    xdft = zeros(N_b, n_mic);
    for M = 1:n_mic
        % fft with hanning window to all microphones, index according to 50% overlap
        xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx = psdx + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx = (8/3)*psdx/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx(2:end-1,:) = 2*psdx(2:end-1,:); % Single-side
fprintf([reverseStr, '\tAll blocks evaluated!\n']);
[sl4, f4] = PSDToSpecLev(mean(psdx,2), fband, Fs);
clear data psdx xdft;

load(datafile5);
data = data(:,mic_sel);
[N, n_mic] = size(data);
t = N/50e3; % sec
N_b = t_b*Fs;
n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
fprintf('\tUsing %d blocks\n', n_blocks);
df = Fs/N_b;
psdx = zeros(floor(N_b/2), n_mic);
freq = (0:floor(N_b/2)-1)*df;
reverseStr = [];
for B = 1:n_blocks
    msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    xdft = zeros(N_b, n_mic);
    for M = 1:n_mic
        % fft with hanning window to all microphones, index according to 50% overlap
        xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx = psdx + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx = (8/3)*psdx/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx(2:end-1,:) = 2*psdx(2:end-1,:); % Single-side
fprintf([reverseStr, '\tAll blocks evaluated!\n']);
[sl5, f5] = PSDToSpecLev(mean(psdx,2), fband, Fs);
clear data psdx xdft;

load(datafile6);
data = data(:,mic_sel);
[N, n_mic] = size(data);
t = N/50e3; % sec
N_b = t_b*Fs;
n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
fprintf('\tUsing %d blocks\n', n_blocks);
df = Fs/N_b;
psdx = zeros(floor(N_b/2), n_mic);
freq = (0:floor(N_b/2)-1)*df;
reverseStr = [];
for B = 1:n_blocks
    msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    xdft = zeros(N_b, n_mic);
    for M = 1:n_mic
        % fft with hanning window to all microphones, index according to 50% overlap
        xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx = psdx + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx = (8/3)*psdx/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx(2:end-1,:) = 2*psdx(2:end-1,:); % Single-side
fprintf([reverseStr, '\tAll blocks evaluated!\n']);
[sl6, f6] = PSDToSpecLev(mean(psdx,2), fband, Fs);
clear data psdx xdft;

figure;
hold on
h1 = plot(f1/1000, sl1, 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
h4 = plot(f4/1000, sl4, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
h7 = plot(fb/1000, slb, 'LineWidth', 1, 'Color', [0.4660, 0.6740, 0.1880]);
hl = legend('No flow', 'Flow', 'BG', 'Location', 'NorthEast');
hold off
plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', 'RPM = 4400 [min$^{-1}$]', [-0.1 5], [0 80], ...
    0:1:5e0, 0:10:80, 20, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\SPECTRUMLEVEL_ARRAY_PROPRPM1']);

figure;
hold on
h2 = plot(f2/1000, sl2, 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
h5 = plot(f5/1000, sl5, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
h7 = plot(fb/1000, slb, 'LineWidth', 1, 'Color', [0.4660, 0.6740, 0.1880]);
hl = legend('No flow', 'Flow', 'BG', 'Location', 'NorthEast');
hold off
plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', 'RPM = 7000 [min$^{-1}$]', [-0.1 5], [0 80], ...
    0:1:5e0, 0:10:80, 20, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\SPECTRUMLEVEL_ARRAY_PROPRPM2']);

figure;
hold on
h3 = plot(f3/1000, sl3, 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
h6 = plot(f6/1000, sl6, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
h7 = plot(fb/1000, slb, 'LineWidth', 1, 'Color', [0.4660, 0.6740, 0.1880]);
hl = legend('No flow', 'Flow', 'BG', 'Location', 'NorthEast');
hold off

plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', 'RPM = 7600 [min$^{-1}$]', [-0.1 5], [0 80], ...
    0:1:5e0, 0:10:80, 20, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\SPECTRUMLEVEL_ARRAY_PROPRPM3']);
if save_data; close; close; close; end


end

function [SpecLev, fi] = PSDToSpecLev(psdx, fband, Fs)
% The frequency 'spacing' is the amount of points times the frequency
% resolution of the psd! After performing apply logarithmic sum of levels
% to obtain SPL (or energy), this in contrast to integration of the PSD. 
% In this code there is no check if points or nbands is not exact integer!

N = 2*length(psdx);
points = fband*N/Fs;

nbands = N/2/points;

SpecLev = zeros(nbands,1);
for I = 1:nbands
    SpecLev(I) = 20*log10(sqrt(sum(psdx((I-1)*points+1:points*I)*Fs/N))/2e-5);
end

fi = (0:fband:fband*(nbands-1)).'+fband/2;

end