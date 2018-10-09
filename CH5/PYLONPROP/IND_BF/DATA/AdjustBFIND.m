function AdjustBFIND
close all;
filenaam = 'PYLPROPBLOW_900-8000.mat';
load(filenaam);
I = 3;
addpath('O:\MATLAB Signal Processing Files');

fileender = {'PROP'; 'PYLPROP'; 'PYLPROPBLOW'};
titles = {'Isolated propeller';
          'Pylon propeller';
          'Pylon propeller blowing'};
strtit = titles{I};
save_data = 1;


%% CB action

dynamic_range = 6;
maxvalxz = max(SPLxzCB(:));

figure;
contourf(x, z, SPLxzCB, (round(maxvalxz)-dynamic_range):0.5/2:round(maxvalxz));

if strcmp(fileender{I}(1:4), 'PROP')
    layoutPylonprop(gca,0);
else
    layoutPylonprop(gca,1);
end

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [strtit], ...
    [xmin xmax], [zmin zmax], linspace(xmin, xmax, 7), ...
    linspace(zmin, zmax, 6), 16, 'on', 'on', 1, ...
    [1 round(maxvalxz)-dynamic_range round(maxvalxz)], ...
    '$L_{\mathrm{p}}$ [dB]', save_data, ...
    ['..\..\BFCBxz' filenaam(1:end-4) '_y=' num2str(scan_plane_Y)]);
% close;

%% CLEAN SC

dynamic_range = 12;
maxvalxz = max(SPLxzSC(:));

figure;
hIm = imagesc(x, z, SPLxzSC, [round(maxvalxz)-dynamic_range round(maxvalxz)]);
hIm.AlphaData = SPLxzSC>(round(maxvalxz)-dynamic_range);

if strcmp(fileender{I}(1:4), 'PROP')
    layoutPylonprop(gca,0);
else
    layoutPylonprop(gca,1);
end

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [strtit], ...
    [xmin xmax], [zmin zmax], linspace(xmin, xmax, 7), ...
    linspace(zmin, zmax, 6), 16, 'on', 'on', 1, ...
    [1 round(maxvalxz)-dynamic_range round(maxvalxz)], ...
    '$L_{\mathrm{p}}$ [dB]', save_data, ...
    ['..\..\BFSCxz' filenaam(1:end-4) '_y=' num2str(scan_plane_Y)]);
% close;

%% HR CLEAN SC

dynamic_range = 12;

Xr = ones(numel(z),1)*x - mean(mic_config(:,1),1);
Zr = z.'*ones(1,numel(x)) - mean(mic_config(:,3),1);
Yr = -1 - mean(mic_config(:,2),1);
Rsep = sqrt(Xr.^2 + Yr.^2 + Zr.^2);

SPLxzHR = SPLxzHR + 20*log10(4*pi*8.35) - 20*log10(4*pi*Rsep);

maxvalxz = max(SPLxzHR(:));

figure;
hIm = imagesc(x, z, SPLxzHR, [round(maxvalxz)-dynamic_range round(maxvalxz)]);
hIm.AlphaData = SPLxzHR>(round(maxvalxz)-dynamic_range);

if strcmp(fileender{I}(1:4), 'PROP')
    layoutPylonprop(gca,0);
else
    layoutPylonprop(gca,1);
end

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [strtit], ...
    [xmin xmax], [zmin zmax], linspace(xmin, xmax, 7), ...
    linspace(zmin, zmax, 6), 16, 'on', 'on', 1, ...
    [1 round(maxvalxz)-dynamic_range round(maxvalxz)], ...
    '$L_{\mathrm{p}}$ [dB]', save_data, ...
    ['..\..\BFHRxz' filenaam(1:end-4) '_y=' num2str(scan_plane_Y)]);
% close;

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
