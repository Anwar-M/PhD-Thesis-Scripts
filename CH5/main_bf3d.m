function main_bf3d
save_data = 0;
pref = 2e-5;
fs = 50e3;

xsource = [0 -1 0];

addpath('O:\MATLAB Signal Processing Files');

main_path = 'O:\APIAN\APIAN Pylon Processed\Calib\';
save_path = 'O:\PhD Thesis\RESULTS\CH5\LOUDSPEAKER\';
remove_mic = [39];

file_path{1} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_09-16-46_nowind';
save_path = [];


[data, mic_config] = read_data(file_path{1}, main_path, remove_mic, 0);
fprintf('\tMean SPL over array: %.2f [dB]\n',20*log10(mean(rms(data,1))/2e-5));
N_mics = size(mic_config, 1);

%% CSM
c = 344.5;
t_start = 0;
t_end = 10;

f_low = 1e3/sqrt(2);
f_high = 2e3*sqrt(2);

% overlap 50%, blocks of 1 sec start at 0 end at 10
[CSM, freqs] = developCSM(data, f_low, f_high, 50e3, 0.01, .5, t_start, t_end);
n_freqs = length(freqs);

%% Beamforming
xmin = -1.5; xmax = 1.5;
ymin = -2; ymax = 1;

reso = 0.05;
rel_range = (-1:reso:2);
z = rel_range+xsource(1,3);
% z = [1:0.1:1.6 xsource(3,3) xsource(2,3) xsource(1,3) 2.1 2.2 3];

tic;
A = [];
for I = 1:length(z)
    msg = sprintf('\tBeamforming plane %d/%d, z = %.2f [m]...\n', I, length(z), z(I));
    fprintf(msg);
    [x, y, B] = FastBeamforming1(CSM, z(I), freqs, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
%     [x, y, B] = FastBeamformingPS(CSM, z(I), freqs, ...
%                 [xmin xmax ymin ymax], reso, mic_config.', c);
%     [x, y, B] = FastBeamforming4(CSM, z(I), freqs, ...
%                 [xmin xmax ymin ymax], reso, mic_config.', c);
%     [x, y, B] = CleanSC(CSM, z(I), freqs, [xmin xmax ymin ymax], reso, ...
%                         mic_config.', c);
    A = cat(3,A,B);
end
disp(toc);

set(0,'defaulttextInterpreter','latex') 
zs = [xsource(1,3)];
xs = [xsource(1,1)];
ys = [xsource(1,2)];

bf3d_script(A, x, y, z, zs, xs, ys, 12, mic_config, 1, 1, xsource);
% title(['$f =$ ' num2str(frequency) ' [Hz]'], 'Interpreter', 'Latex');
set(gca, 'FontSize', 14);
set(gca,'TickLabelInterpreter','latex');
keyboard;

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
