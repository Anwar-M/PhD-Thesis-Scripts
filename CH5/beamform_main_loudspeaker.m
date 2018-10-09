function beamform_main_loudspeaker
addpath('O:\MATLAB Signal Processing Files');
main_path = 'O:\APIAN\APIAN Pylon Processed\Calib\';
save_path = 'O:\PhD Thesis\RESULTS\CH5\LOUDSPEAKER\';
save_data = 0;
XS = [0 -1 0];
remove_mic = [39];
freq_lb = [1000/sqrt(2) 1000/sqrt(2) 4000/sqrt(2) 4000/sqrt(2) 1000/sqrt(2) 4000/sqrt(2)];
freq_ub = [2000*sqrt(2) 2000*sqrt(2) 8000*sqrt(2) 8000*sqrt(2) 2000*sqrt(2) 8000*sqrt(2)];

titles = {'1 and 2 [kHz] octave band';
          '1 and 2 [kHz] octave band';
          '4 and 8 [kHz] octave band';
          '4 and 8 [kHz] octave band';
          '1 and 2 [kHz] octave band';
          '4 and 8 [kHz] octave band'};
fileender = {'noflow'; 'flow'; 'noflow'; 'flow'; 'flow'; 'flow';};

t_start = 0;
t_end = 10;
Tblock = 0.01;

% 1 & 2 kHz octave no wind
file_path{1} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_09-16-46_nowind';
% 1 & 2 kHz octave
file_path{2} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_11-29-19';
% 4 & 8 kHz octave no wind
file_path{3} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_09-28-02_nowind';
% 4 & 8 kHz octave
file_path{4} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_11-59-10';


% 1 & 2 kHz octave wind 40
file_path{5} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_09-51-44_U40';
% 4 & 8 kHz octave wind 40
file_path{6} = 'O:\APIAN\APIAN Pylon Processed\Calib\2014-09-29_10-01-53_U40';

for I = [3];%1:4
    [data, mic_pos] = read_data(file_path{I}, main_path, remove_mic, 0);
    fprintf('\tMean SPL over array: %.2f [dB]\n',20*log10(mean(rms(data,1))/2e-5));
    beamform(data, mic_pos, freq_lb(I), freq_ub(I), Tblock, t_start, t_end, XS, titles{I}, save_data, save_path, fileender{I});
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

function beamform(data, mic_config, flb, fub, Tblock, t_start, t_end, XS, strtit, save_data, save_path, fileender)
Fs = 50e3;
% Sound speed
c = 344.5;

reso = 0.01;
xmin = -1.5; xmax = 1.5;
zmin = -1.5; zmax = 1.5;

ymin = -2.5; ymax = 0.5;

% overlap 50%, blocks of 1 sec start at 0 end at 10
[CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);

scan_plane_Y = [-1.0];
scan_plane_Z = 0;

% Beamform, Sarradj form 3. Get values at array center
[x, z, Axz] = FastBeamforming3(CSM, scan_plane_Y , freqs, ...
        [xmin xmax zmin zmax], reso, [mic_config(:,1).';mic_config(:,3).';mic_config(:,2).'], c);
[x, y, Axy] = FastBeamforming3(CSM, scan_plane_Z , freqs, ...
        [xmin xmax ymin ymax], reso, mic_config.', c);
% [x, z, Axz] = FastBeamforming3Conv(CSM, scan_plane_Y , freqs, ...
%         [xmin xmax zmin zmax], reso, [mic_config(:,1).';mic_config(:,3).';mic_config(:,2).'], c, [40 0 0]*(8.35-4.6)/8.35);
% [x, y, Axy] = FastBeamforming3Conv(CSM, scan_plane_Z , freqs, ...
%         [xmin xmax ymin ymax], reso, mic_config.', c, [40 0 0]*(8.35-4.6)/8.35);
SPLxz = 20*log10( sqrt(real(Axz)) / 2e-5 );
SPLxy = 20*log10( sqrt(real(Axy)) / 2e-5 );

dynamic_range = 6;
maxvalxz = max(SPLxz(:));

figure;
contourf(x, z, SPLxz, (round(maxvalxz)-dynamic_range):0.5:round(maxvalxz));
hold on
h1 = plot(XS(1),XS(3),'kx');
set(h1,'MarkerSize',10)
set(h1,'LineWidth',2)
hold off

[dummy, vert_ind] = max(SPLxz);
[~, hor_ind] = max(dummy);
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
Xsh = xlim(1) + (hor_ind-1)*reso;
Zsh = ylim(1) + (vert_ind(hor_ind)-1)*reso;
placement_x = 1*(xlim(2)-xlim(1))/40 + xlim(1);
placement_y = 1*(ylim(2)-ylim(1))/30 + ylim(1);

line([xlim(1); xlim(2)], [Zsh; Zsh], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([Xsh; Xsh], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);

text(placement_x, placement_y, ['Max: (' num2str(Xsh) ', ' num2str(Zsh) ')'], ...
    'Color','k','FontSize', 14,'Interpreter','LaTex', 'BackgroundColor', 0.93*[1 1 1]);

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [strtit], ...
    [xmin xmax], [zmin zmax], linspace(xmin, xmax, 5), linspace(zmin, zmax, 5), 16, ...
    'on', 'on', 1, [1 round(maxvalxz)-dynamic_range round(maxvalxz)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path 'BFxz' fileender '_y=' num2str(scan_plane_Y) '_fl=' num2str(round(flb))]);

maxvalxy = max(SPLxy(:));
figure;
contourf(x, y, SPLxy, (round(maxvalxy)-dynamic_range):0.5:round(maxvalxy));
hold on
h1 = plot(XS(1),XS(2),'kx');
set(h1,'MarkerSize',10)
set(h1,'LineWidth',2)
hold off

[dummy, vert_ind] = max(SPLxy);
[~, hor_ind] = max(dummy);
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
Xsh = xlim(1) + (hor_ind-1)*reso;
Ysh = ylim(1) + (vert_ind(hor_ind)-1)*reso;
placement_x = 1*(xlim(2)-xlim(1))/40 + xlim(1);
placement_y = 1*(ylim(2)-ylim(1))/30 + ylim(1);

line([xlim(1); xlim(2)], [Ysh; Ysh], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([Xsh; Xsh], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);

text(placement_x, placement_y, ['Max: (' num2str(Xsh) ', ' num2str(Ysh) ')'], ...
    'Color','k','FontSize', 14,'Interpreter','LaTex', 'BackgroundColor', 0.93*[1 1 1]);
 
plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [strtit], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 round(maxvalxy)-dynamic_range round(maxvalxy)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path 'BFxy' fileender '_z=' num2str(scan_plane_Z) '_' num2str(round(flb))]);



