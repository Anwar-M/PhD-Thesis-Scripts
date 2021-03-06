function main_spectrumlevel_loudspeaker
close all;
save_data = 1;

main_path = 'O:\APIAN\APIAN Pylon Processed\Calib\';
addpath('O:\MATLAB Signal Processing Files');
mainresultpath = 'O:\PhD Thesis\RESULTS\CH5\LOUDSPEAKER';

dirdata{1} = '2014-09-29_09-16-46_nowind';
dirdata{2} = '2014-09-29_09-28-02_nowind';
dirdata{3} = '2014-09-29_11-29-19';
dirdata{4} = '2014-09-29_11-59-10';

% mic_sel = 1:8;
Fs = 50e3; % sample frequency
t_b = 0.1;
fband = 50;

for I = 1:4
    file_path = [main_path dirdata{I}];
    
    remove_mic = [39];
    [data, ~] = read_data(file_path, main_path, remove_mic, 0);
%     data = data(:,mic_sel);

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
    
    [sl(:,I), f1(:,I)] = PSDToSpecLev(mean(psdx1,2), fband, Fs);

    clear data psdx1 xdft;

end
hold on
plot(f1(:,1)/1000, sl(:,1), 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
plot(f1(:,3)/1000, sl(:,3), 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
line([1;1]./sqrt(2),[0;100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
line(2*sqrt(2)*[1;1],[0;100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hl = legend('No flow', '$U$ = 60 [m/s]', 'Location', 'NorthEast');
set(hl, 'Interpreter', 'LaTex');
plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', ['1 and 2 [kHz] octave band'], [-0.5 25], [0 100], ...
    0:5:25e0, 0:20:100, 20, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\SPECTRUMLEVEL_LS_1_2_KHZ']);
hold off

figure
hold on
plot(f1(:,2)/1000, sl(:,2), 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
plot(f1(:,4)/1000, sl(:,4), 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
line([4;4]./sqrt(2),[0;100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
line(8*sqrt(2)*[1;1],[0;100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hl = legend('No flow', '$U$ = 60 [m/s]', 'Location', 'NorthEast');
set(hl, 'Interpreter', 'LaTex');
plot_settings_font(gca, '$f$ [kHz]', 'Spectrum Level [dB]', ['4 and 8 [kHz] octave band'], [-0.5 25], [0 100], ...
    0:5:25e0, 0:20:100, 20, 'on', 'on', 0, 0, [], save_data, ...
    [mainresultpath '\SPECTRUMLEVEL_LS_4_8_KHZ']);
hold off

end

function [data, ARRAY_ROT] = read_data(file_path, main_path, remove_mic, calib)
fprintf('Read in measurement data\n')

if calib
    fprintf('Loading LabVIEW calibrated data...\n');
    fid = fopen([file_path '\acoustic_data.cal'],'r');
    data = fread(fid, [64 inf], 'single', 0, 'l').';
    fclose(fid);
    fprintf('LabVIEW calibrated data loaded!\n');
else
    fprintf('Loading raw data...\n');
    fid = fopen([file_path '\acoustic_data'],'r');
    data = fread(fid,[64 inf],'int16',0,'l').';
    fclose(fid);
    fprintf('Performing linear MATLAB calibration...\n');
    data = data*2.5/32768/3.985;
    calib_file = 'O:\APIAN\APIAN Pylon Processed\TUDArrayCalib\Microphone calibration\miccalibBundle.txt';
    mic_responses = importdata(calib_file,'\t',0);
    mic_responses(:, [1 3 4]) = [];

    array4mic_order = 8*(1:64) + 56 - 63 * ceil((1:64)/8);
    data(:, array4mic_order) = data(:, array4mic_order)./(ones(size(data,1),1)*mic_responses(1:64).');
    data = data - (ones(size(data,1),1)*mean(data));
    fprintf('Raw data calibrated!\n');
end

load([main_path '..\channels.mat']);
load([main_path '..\ARRAY_ROT_fix.mat']);
ARRAY_ROT = ARRAY_ROT(channels(:,1), :);

if ~isempty(remove_mic)
    data(:, remove_mic) = [];
    ARRAY_ROT(channels(remove_mic,1), :) = [];
end

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