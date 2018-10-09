function beamform_main_pylonprop
% close all;
save_data = 1;

main_path = 'O:\APIAN\APIAN Pylon Processed\Pylon-Prop\';
addpath('O:\MATLAB Signal Processing Files');
mainresultpath = 'O:\PhD Thesis\RESULTS\CH5\PYLONPROP';

dirdata{1} = 'IsolatedProp';
dirdata{2} = 'PylonPropNoBlow';
dirdata{3} = 'PylonPropBlow';

remove_mic = [39];
freq_lb = 900*[1 1 1];
freq_ub = 8000*[1 1 1];

titles = {'Isolated propeller';
          'Pylon propeller';
          'Pylon propeller blowing'};
fileender = {'PROP'; 'PYLPROP'; 'PYLPROPBLOW'};

t_start = 0;
t_end = 32;
Tblock = 0.25;

for I = 1:3
    file_path = [main_path dirdata{I}];
    [data, mic_pos] = read_data(file_path, main_path, remove_mic, 0);
    fprintf('\tMean SPL over array: %.2f [dB]\n',20*log10(mean(rms(data,1))/2e-5));
    beamform(data, mic_pos, freq_lb(I), freq_ub(I), Tblock, t_start, t_end, titles{I}, save_data, mainresultpath, fileender{I});
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

% remove diag
n_mic = size(CSM,1);
n_f = size(CSM,3);
CSMvec = CSM(:);
for ff = 1:n_f
    CSMvec((n_mic+1)*(0:(n_mic-1))+1+n_mic*n_mic*(ff-1)) = 0;
end
CSM = reshape(CSMvec,n_mic,n_mic,n_f);
clear CSMvec;

scan_plane_Y = -1.0;
scan_plane_Z = -0.97;

% Beamform, Sarradj form 3. Get values at array center
% [x, z, Axz] = FastBeamforming3(CSM, scan_plane_Y , freqs, ...
%         [xmin xmax zmin zmax], reso, [mic_config(:,1).';mic_config(:,3).';mic_config(:,2).'], c);
% [x, y, Axy] = FastBeamforming3(CSM, scan_plane_Z , freqs, ...
%         [xmin xmax ymin ymax], reso, mic_config.', c);
[x, z, Axz] = FastBeamforming3Conv(CSM, scan_plane_Y , freqs, ...
        [xmin xmax zmin zmax], reso, [mic_config(:,1).';mic_config(:,3).';mic_config(:,2).'], c, [60 0 0]*(8.35-4.6)/8.35);
% [x, y, Axy] = FastBeamforming3Conv(CSM, scan_plane_Z , freqs, ...
%         [xmin xmax ymin ymax], reso, mic_config.', c, [60 0 0]*(8.35-4.6)/8.35);

Axz((real(Axz)<0)) = 0;
% Axy((real(Axy)<0)) = 0;
    
SPLxz = 20*log10( sqrt(real(Axz)) / 2e-5 );
% SPLxy = 20*log10( sqrt(real(Axy)) / 2e-5 );

dynamic_range = 6;
maxvalxz = max(SPLxz(:));

if save_data
    file_data = [fileender 'dr.mat'];
    save([save_path '\DATA\' file_data], 'Fs', 'c', 'reso', 'xmin', 'xmax', 'zmin', 'zmax', ...
        'CSM', 'freqs', 'Axz', 'scan_plane_Y', 'mic_config', 'flb', 'fub', 'Tblock', 't_start', 't_end', ...
        'x', 'z', 'SPLxz', 'dynamic_range', 'maxvalxz');
end

figure;
contourf(x, z, SPLxz, (round(maxvalxz)-dynamic_range):0.5:round(maxvalxz));
layoutPylonprop(gca, 1);

% [dummy, vert_ind] = max(SPLxz);
% [~, hor_ind] = max(dummy);
% xlim = get(gca, 'XLim');
% ylim = get(gca, 'YLim');
% Xsh = xlim(1) + (hor_ind-1)*reso;
% Zsh = ylim(1) + (vert_ind(hor_ind)-1)*reso;
% placement_x = 1*(xlim(2)-xlim(1))/40 + xlim(1);
% placement_y = 1*(ylim(2)-ylim(1))/30 + ylim(1);
% 
% line([xlim(1); xlim(2)], [Zsh; Zsh], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
% line([Xsh; Xsh], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
% 
% text(placement_x, placement_y, ['Max: (' num2str(Xsh) ', ' num2str(Zsh) ')'], ...
%     'Color','k','FontSize', 14,'Interpreter','LaTex', 'BackgroundColor', 0.93*[1 1 1]);

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [strtit], ...
    [xmin xmax], [zmin zmax], linspace(xmin, xmax, 5), linspace(zmin, zmax, 5), 16, ...
    'on', 'on', 1, [1 round(maxvalxz)-dynamic_range round(maxvalxz)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\BFxz' fileender 'dr_y=' num2str(scan_plane_Y)]);

% maxvalxy = max(SPLxy(:));
% figure;
% contourf(x, y, SPLxy, (round(maxvalxy)-dynamic_range):0.5:round(maxvalxy));
% 
% layoutPylonpropxy(gca,1);
% 
% [dummy, vert_ind] = max(SPLxy);
% [~, hor_ind] = max(dummy);
% xlim = get(gca, 'XLim');
% ylim = get(gca, 'YLim');
% Xsh = xlim(1) + (hor_ind-1)*reso;
% Ysh = ylim(1) + (vert_ind(hor_ind)-1)*reso;
% placement_x = 1*(xlim(2)-xlim(1))/40 + xlim(1);
% placement_y = 1*(ylim(2)-ylim(1))/30 + ylim(1);
% 
% line([xlim(1); xlim(2)], [Ysh; Ysh], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
% line([Xsh; Xsh], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
% 
% text(placement_x, placement_y, ['Max: (' num2str(Xsh) ', ' num2str(Ysh) ')'], ...
%     'Color','k','FontSize', 14,'Interpreter','LaTex', 'BackgroundColor', 0.93*[1 1 1]);
%  
% plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [strtit], ...
%     [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
%     'on', 'on', 1, [1 round(maxvalxy)-dynamic_range round(maxvalxy)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\BFCBxy' fileender '_z=' num2str(round(scan_plane_Z)) '_' num2str(round(flb))]);

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

function layoutPylonpropxy(hAxes, pylon)
% edits @ 24-9-2018

xp = 0:1/100:97/100;
yp = 5*0.1*(0.2969*sqrt(xp)-0.126*xp-0.3516*xp.^2+0.2843*xp.^3-0.1015*xp.^4);
Lp = 0.489;

% Pylon-prop config
parkerpen = 8;

% Prop
Xp = [-.37 -1.0 -.4];
Xpu = [-.37 -1.2540 -0.4];
Xpd = [-.37 -0.7460 -0.4];

% Pylon
Xpyl_corn1 = [-0.5820 -1.0 -.4];
Xpyl_corn2 = [-0.5820 -1.0 -1.3];
Xpyl_corn3 = [-1.0710 -1.0 -1.3];
Xpyl_corn4 = [-1.0710 -1.0 -.4];

% DNW support above
Sup1 = [Xp + [.096 0 0.09]; [-0.34 -1 -0.3425]; [-0.34 -1 -0.4575]; Xp + [.096 0 -0.09]; ...
        Xp + [1.461 0 -0.09]; 
        [1.395 -1 -0.21-.383-.017]; [1.796 -1 -0.25-.383-.017]; [1.796 -1 0.25-.383-.017]; [1.395 -1 0.21-.383-.017];...
        Xp + [1.461 0 0.09]];
    
% DNW support below
Sup2 = [Xpyl_corn3; Xpyl_corn3 + [1.167 0 0]; Xpyl_corn3 + [1.167+1.7 0 .25-.083];
        Xpyl_corn3 + [1.167+1.7 0 -.25-.083]; Xpyl_corn3 + [1.167 0 -.166]; ...
        Xpyl_corn3 + [-.333 0 -.166]; Xpyl_corn3 + [-.333-.075 0 -.12]; ...
        Xpyl_corn3 + [-.433 0 -.083]; 
        Xpyl_corn3 + [-.333-.075 0 -.04];...
        Xpyl_corn3 + [-.333 0 0];];

hold(hAxes, 'on')
% plot(hAxes, Xp(1),Xp(2),'ko','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpu(1),Xpu(2),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpd(1),Xpd(2),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpd(1),Xpd(2),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);

plot(hAxes, Lp*xp+Xpyl_corn3(1), Lp*yp+Xpyl_corn3(2), 'k', Lp*xp+Xpyl_corn3(1), -Lp*yp+Xpyl_corn3(2), 'k', 'LineWidth', 1.5);

xpr = [-0.330000000000000;-0.420000000000000;-0.430000000000000;-0.420000000000000;-0.505000000000000-0.02;-0.550000000000000-0.02;-0.550000000000000-0.02;-0.505000000000000-0.02;-0.420000000000000;-0.430000000000000;-0.420000000000000;-0.330000000000000;-0.350000000000000;-0.350000000000000];
zpr = [-0.146000000000000;-0.245000000000000;-0.300000000000000;-0.345000000000000;-0.360000000000000;-0.385000000000000;-0.415000000000000;-0.440000000000000;-0.455000000000000;-0.500000000000000;-0.555000000000000;-0.654000000000000;-0.470000000000000;-0.330000000000000];

patch(xpr+0.01,zpr-0.6, 'w', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
% patch([xpr(4)+0.01; xpr(4)+0.01],[zpr(4) zpr(9)], 'w', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% if pylon
% patch([Xpyl_corn1(1); Xpyl_corn2(1); Xpyl_corn3(1); Xpyl_corn4(1); Xpyl_corn1(1)], ...
%       [Xpyl_corn1(3); Xpyl_corn2(3); Xpyl_corn3(3); Xpyl_corn4(3); Xpyl_corn1(3)], 'k', ...
%       'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
% end
% 
patch(Sup1(:,1), Sup1(:,3)-0.6, 'k--', 'Parent', hAxes, 'FaceColor', 'none', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
patch(Sup2(:,1), Sup2(:,3)+.383, 'k--', 'Parent', hAxes, 'FaceColor', 'none', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
hold(hAxes, 'off')

