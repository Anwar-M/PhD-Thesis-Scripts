function beamformVar_main_pylonprop_indfreqHRFIX
% There was an error with the steering vector for adap hr clean sc. Is now
% fixed and maps for only hr needed to be redone. Immediate correction for
% spl required at array center.
save_data = 1;

main_path = 'O:\APIAN\APIAN Pylon Processed\Pylon-Prop\';
addpath('O:\MATLAB Signal Processing Files');
mainresultpath = 'O:\PhD Thesis\RESULTS\CH5\PYLONPROP\IND_BF';

dirdata{1} = 'IsolatedProp';
dirdata{2} = 'PylonPropNoBlow';
dirdata{3} = 'PylonPropBlow';

remove_mic = [39];
freq_lb = [475 995 1500 2005 2520 3025 3535 4040 1332 1610 900];
freq_ub = [525 1045 1550 2055 2570 3075 3585 4090 1382 1900 8000];

titles = {'Isolated propeller';
          'Pylon propeller';
          'Pylon propeller blowing'};
fileender = {'PROP'; 'PYLPROP'; 'PYLPROPBLOW'};

t_start = 0;
t_end = 32;
Tblock = 0.1;

for I = 1:3
    file_path = [main_path dirdata{I}];
    [data, mic_pos] = read_data(file_path, main_path, remove_mic, 0);
    fprintf('\tMean SPL over array: %.2f [dB]\n',20*log10(mean(rms(data,1))/2e-5));
    for F = 1:(numel(freq_lb)-1)
        fprintf('\t%d - %d [Hz]\n', freq_lb(F), freq_ub(F));
        beamform(data, mic_pos, freq_lb(F), freq_ub(F), Tblock, t_start, ...
            t_end, titles{I}, save_data, mainresultpath, ...
            [fileender{I} '_' num2str(freq_lb(F)) '-' num2str(freq_ub(F))]);
    end
end

fprintf('\tStart with all freqs (900-8000 Hz)...\n');

t_start = 0;
t_end = 20;
Tblock = 0.025;
for I = 1:3
    file_path = [main_path dirdata{I}];
    [data, mic_pos] = read_data(file_path, main_path, remove_mic, 0);
    F = 11;
    fprintf('\t%d - %d [Hz]\n', freq_lb(F), freq_ub(F));
    beamform(data, mic_pos, freq_lb(F), freq_ub(F), Tblock, t_start, ...
        t_end, titles{I}, save_data, mainresultpath, ...
        [fileender{I} '_' num2str(freq_lb(F)) '-' num2str(freq_ub(F))]);
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

function beamform(data, mic_config, flb, fub, Tblock, t_start, t_end, strtit, save_data, save_path, fileender)
Fs = 50e3;
% Sound speed
c = 344.5;

reso = 0.01;
xmin = -1.5; xmax = 1.5;
zmin = -2; zmax = 0.5;

ymin = -2.5; ymax = 0.5;

% overlap 50%, blocks of 1 sec start at 0 end at 10
[CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);

scan_plane_Y = -1.0;
    
[x, z, AxzHR, ~] = adaptive_HR_CleanSCConv_mod(CSM, scan_plane_Y , freqs, ...
         [xmin xmax zmin zmax], reso, ...
         [mic_config(:,1).';mic_config(:,3).';mic_config(:,2).'], c, ...
         [60 0 0]*(8.35-4.6)/8.35, 2);
     
Xr = ones(numel(z),1)*x - mean(mic_config(:,1),1);
Zr = z.'*ones(1,numel(x)) - mean(mic_config(:,3),1);
Yr = -1 - mean(mic_config(:,2),1);
Rsep = sqrt(Xr.^2 + Yr.^2 + Zr.^2);
    
SPLxzHR = 20*log10( sqrt(real(AxzHR)) ./ (4*pi*2e-5*Rsep));

if save_data
    file_data = [fileender 'HRFIX.mat'];
    save([save_path '\DATA\' file_data], 'Fs', 'c', 'reso', 'xmin', 'xmax', 'zmin', 'zmax', ...
        'freqs', 'AxzHR', 'scan_plane_Y', 'mic_config', 'flb', 'fub', 'Tblock', 't_start', 't_end', ...
        'x', 'z', 'SPLxzHR', 'Rsep');
end

%% HR CLEAN SC

dynamic_range = 12;
maxvalxz = max(SPLxzHR(:));

figure;
hIm = imagesc(x, z, SPLxzHR, [round(maxvalxz)-dynamic_range round(maxvalxz)]);
hIm.AlphaData = SPLxzHR>(round(maxvalxz)-dynamic_range);

if strcmp(fileender(1:4), 'PROP')
    layoutPylonprop(gca,0);
else
    layoutPylonprop(gca,1);
end

xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
placement_x = 13*(xlim(2)-xlim(1))/24 + xlim(1);
placement_y = 5*(ylim(2)-ylim(1))/120 + ylim(1);
placement_x2 = 43*(xlim(2)-xlim(1))/48 + xlim(1);
placement_y2 = 112*(ylim(2)-ylim(1))/120 + ylim(1);
text(placement_x, placement_y, [num2str(flb) ' - ' num2str(fub) ' [Hz]'], ...
    'Color','k','FontSize', 16,'Interpreter','LaTex', 'BackgroundColor', [1 1 1]);
text(placement_x2, placement_y2, ['HR'], ...
    'Color','k','FontSize', 16,'Interpreter','LaTex', 'BackgroundColor', [1 1 1]);

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [strtit], ...
    [xmin xmax], [zmin zmax], linspace(xmin, xmax, 7), ...
    linspace(zmin, zmax, 6), 16, 'on', 'on', 1, ...
    [1 round(maxvalxz)-dynamic_range round(maxvalxz)], ...
    '$L_{\mathrm{p}}$ [dB]', save_data, ...
    [save_path '\BFHRxz' fileender '_y=' num2str(scan_plane_Y)]);
close;

function layoutPylonprop(hAxes, pylon)
% edits @ 19-9-2018

% Pylon-prop config
parkerpen = 8;

% Prop
Xp = [-.37 -1.0 -.4];
Xpu = [-.37 -1.0 -0.146];
Xpd = [-.37 -1.0 -0.654];

% Pylon
Xpyl_corn1 = [-0.5820 -1.0 -.4];
Xpyl_corn2 = [-0.5820 -1.0 -1.3];
Xpyl_corn3 = [-1.0710 -1.0 -1.3];
Xpyl_corn4 = [-1.0710 -1.0 -.4];

% DNW support above
Sup1 = [Xp + [.096 0 0.09]; [-0.34 -1 -0.3425]; [-0.34 -1 -0.4575]; Xp + [.096 0 -0.09]; ...
        Xp + [1.461 0 -0.09]; Xp + [1.461 0 0.09]; ...
        Xp + [1.422 0 0.09]; Xp + [1.422 0 0.584-0.214]; Xp + [1.656 0 0.584-0.214]; ...
        Xp + [1.656+0.097 0 0.584-0.214-0.097]; [1.14 -1 -0.36]; [1.091 -1 -0.36]; [1.091 -1 -0.44]; ...
        [1.16 -1 -0.44]; Xp + [1.656+0.097+0.06 0 0.584-0.214-0.097-0.05]; [1.796 -1 -1.133]; 
        Xp + [1.133+2 0 0.584-0.214]; Xp + [1.133+2 0 0.584]; ...
        Xp + [1.130-0.078 0 0.584]; Xp + [1.130-0.078 0 0.584-0.214]; Xp + [1.130 0 0.584-0.214]; Xp + [1.130 0 0.09]];
    
% DNW support below
Sup2 = [Xpyl_corn3; Xpyl_corn3 + [1.167 0 0]; Xpyl_corn3 + [1.167+1.7 0 .25-.083];
        Xpyl_corn3 + [1.167+1.7 0 -.25-.083]; Xpyl_corn3 + [1.167 0 -.166]; ...
        Xpyl_corn3 + [-.333 0 -.166]; Xpyl_corn3 + [-.333-.075 0 -.12]; ...
        Xpyl_corn3 + [-.433 0 -.083]; Xpyl_corn3 + [-.333-.075 0 -.04];...
        Xpyl_corn3 + [-.333 0 0];];

hold(hAxes, 'on')
% plot(hAxes, Xp(1),Xp(3),'ko','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpu(1),Xpu(3),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpd(1),Xpd(3),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);

% patch([Xpu(1); Xpd(1)], [Xpu(3); Xpd(3)], 'k', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
% patch([Xpu(1)+0.05; Xpd(1)], [Xpu(3); Xpd(3)], 'w', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

xpr = [-0.330000000000000;-0.420000000000000;-0.430000000000000;-0.420000000000000;-0.505000000000000-0.02;-0.550000000000000-0.02;-0.550000000000000-0.02;-0.505000000000000-0.02;-0.420000000000000;-0.430000000000000;-0.420000000000000;-0.330000000000000;-0.350000000000000;-0.350000000000000];
zpr = [-0.146000000000000;-0.245000000000000;-0.300000000000000;-0.345000000000000;-0.360000000000000;-0.385000000000000;-0.415000000000000;-0.440000000000000;-0.455000000000000;-0.500000000000000;-0.555000000000000;-0.654000000000000;-0.470000000000000;-0.330000000000000];

patch(xpr+0.01,zpr, 'w', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

patch([xpr(4)+0.01; xpr(4)+0.01],[zpr(4) zpr(9)], 'w', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

if pylon
patch([Xpyl_corn1(1); Xpyl_corn2(1); Xpyl_corn3(1); Xpyl_corn4(1); Xpyl_corn1(1)], ...
      [Xpyl_corn1(3); Xpyl_corn2(3); Xpyl_corn3(3); Xpyl_corn4(3); Xpyl_corn1(3)], 'k', ...
      'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
end
patch(Sup1(:,1), Sup1(:,3), 'k--', 'Parent', hAxes, 'FaceColor', 'none', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
patch(Sup2(:,1), Sup2(:,3), 'k--', 'Parent', hAxes, 'FaceColor', 'none', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
hold(hAxes, 'off')
