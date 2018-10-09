function main_spectrumlevel_4to6
save_data = 1;
include_background = 1;
addpath('O:\MATLAB Signal Processing Files');

datafile1 = 'O:\V-Tunnel 30-10 Plate-Source\Data MAT\Mes4.mat';
datafile2 = 'O:\V-Tunnel 30-10 Plate-Source\Data MAT\Mes5.mat';
datafile3 = 'O:\V-Tunnel 30-10 Plate-Source\Data MAT\Mes6.mat';
mainresultpath = 'O:\PhD Thesis\RESULTS\CH6\OMNISOURCE';
backgroundfile = 'O:\V-Tunnel 30-10 Plate-Source\Data MAT\background.mat';

load('O:\V-Tunnel 30-10 Plate-Source\mic_poses_optim.mat');
mic_config = mic_poses.'; clear mic_poses;
mic_config(:,3) = 0.02;

mic_sel = 41:45;
Fs = 50e3; % sample frequency
t_b = 0.1;
fband = 50;

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

hold on
plot(f1/1000, sl1, 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);

plot(f2/1000, sl2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);

plot(f3/1000, sl3, 'LineWidth', 1, 'Color', [0.4940, 0.1840, 0.5560]);

if include_background&&exist(backgroundfile,'file')
    plot(fb/1000, slb, 'LineWidth', 1, 'Color', [0.4660, 0.6740, 0.1880]);
    hl = legend('Source only', '$d_{\mathrm{object}}$ = 0.5 [m]', '$d_{\mathrm{object}}$ = 0.2 [m]', 'Background', 'Location', 'NorthEast');
    set(hl, 'Interpreter', 'LaTex');
else
    hl = legend('Source only', '$d_{\mathrm{object}}$ = 0.5 [m]', '$d_{\mathrm{object}}$ = 0.2 [m]', 'Location', 'NorthEast');
    set(hl, 'Interpreter', 'LaTex');
end

hold off

plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', ['$d_{\mathrm{array}}$ = 2.03 [m]'], [0 25], [0 50], ...
    0:5:25e0, 0:10:60, 20, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\SPECTRUMLEVEL_ARRAY_4TO6']);
if save_data; close; end


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