function CoherencySpectrumLevels
save_data = 1;
include_background = 0;
addpath('O:\MATLAB Signal Processing Files');

file_path = ['O:\V-Tunnel Coherent Study\DATA MAT\combi_seperate.mat'];
load(file_path);

mainresultpath = 'O:\PhD Thesis\RESULTS\CH2\COHERENCY';
backgroundfile = 'O:\V-Tunnel 13-12 Thrust Experiment\results 15-12\MAT DATA FINAL\background.mat';

load('O:\V-Tunnel Coherent Study\mic_poses_optim.mat');

mic_config = mic_poses.';
clear mic_poses;

Fs = 50e3; % sample frequency
t_b = 0.01;
fband = 100;

fprintf(['\tStart PSD calculations for signals...\n']);
clear data psdxb xdft;

dataL(50e3+1:end,:) = [];
[N, n_mic] = size(dataL);
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
        xdft(:,M) = fft(dataL( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx1 = psdx1 + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx1 = (8/3)*psdx1/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx1(2:end-1,:) = 2*psdx1(2:end-1,:); % Single-side

fprintf([reverseStr, '\tAll blocks evaluated!\n']);
    
[sl1, f1] = PSDToSpecLev(mean(psdx1,2), fband, Fs);

clear dataL psdx1 xdft;

dataR(50e3+1:end,:) = [];
[N, n_mic] = size(dataR);
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
        xdft(:,M) = fft(dataR( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
            hann(N_b));
    end
    psdx2 = psdx2 + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

end
psdx2 = (8/3)*psdx2/n_blocks; % Hanning AMPLITUDE CORRECTION factor
psdx2(2:end-1,:) = 2*psdx2(2:end-1,:); % Single-side

fprintf([reverseStr, '\tAll blocks evaluated!\n']);

[sl2, f2] = PSDToSpecLev(mean(psdx2,2), fband, Fs);

clear dataR psdx2 xdft;

figure;
hold on
plot(f1/1000, sl1, 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);

plot(f2/1000, sl2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
% if include_background&&exist(backgroundfile,'file')
%     plot(freqb/1000, mean_logpsdb, 'LineWidth', 1, 'Color', [0.4660, 0.6740, 0.1880]);
    hl = legend('Left speaker', 'Right speaker', 'Location', 'NorthEast');
    set(hl, 'Interpreter', 'LaTex');
% else
%     hl = legend('Duct on', 'Duct off', 'Location', 'NorthEast');
%     set(hl, 'Interpreter', 'LaTex');
% end

hold off

plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', [], [0 7], [0 40], ...
    0:1e0:7e0, 0:5:40, 16, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\AVG_PSD_ARRAY']);
% close;


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