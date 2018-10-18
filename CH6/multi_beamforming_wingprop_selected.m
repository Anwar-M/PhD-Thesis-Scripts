function multi_beamforming_wingprop_selected
% close all;
save_data = 1;
remdiag = 1;

source_dist = 1.46;
source_height = 1.56;

wing = 1;
wing_dist = 1.071;
wing_heightTE = 1.72;

xs = [0; source_height-1.572; source_dist];

Fs = 50e3;
% Sound speed
c = 345;

addpath('O:\MATLAB Signal Processing Files');
file_path1 = ['O:\V-Tunnel 12-12 Wing-Prop\Data MAT\Meas6.mat'];
file_path2 = ['O:\V-Tunnel 12-12 Wing-Prop\Data MAT\Meas5.mat'];
file_path3 = ['O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\Meas6.mat'];
file_path4 = ['O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\Meas5.mat'];
save_path = [];

load(file_path1);
load('O:\V-Tunnel 12-12 Wing-Prop\mic_poses_optim.mat');
fprintf('\tMeasurement time: %.1f s\n', size(data,1)/Fs);
fprintf('\tAverage SPL over array: %.1f dB\n', 20*log10(mean(rms(data,1))/2e-5));

mic_config = mic_poses.'; clear mic_poses;
mic_config(:,3) = 0.02;
N_mics = size(mic_config, 1);

addpath('O:\MATLAB Signal Processing Files');
mainresultpath = 'O:\PhD Thesis\RESULTS\CH6\PROPWING';

flb = 1275;
fub = 1325;

t_start = 0;
t_end = 30;
Tblock = 0.5;

reso = 0.01;
xmin = -.85; xmax = .85;
ymin = -1; ymax = 0.7;
% z = (~wing)*xs(3) + wing*wing_dist;
% 
% % overlap 50%, blocks of 1 sec start at 0 end at 10
% [CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);
% if remdiag
%     CSM = removeDiagonal(CSM);
% end
% 
% % Beamform, Sarradj form 3. Get values at array center
% [x, y, A] = FastBeamforming3(CSM, z, freqs, ...
%                 [xmin xmax ymin ymax], reso, mic_config.', c);
% % [x, y, Axy] = FastBeamforming3Conv(CSM, scan_plane_Z , freqs, ...
% %         [xmin xmax ymin ymax], reso, mic_config.', c, [60 0 0]*(8.35-4.6)/8.35);
% 
% if remdiag
%     A((real(A)<0)) = 1e-20;
% end
%     
% SPL = 20*log10( sqrt(real(A)) / 2e-5 );
% 
% dynamic_range = 12;
% maxval = max(SPL(:));
% 
% % if save_data
% %     file_data = [fileender 'dr.mat'];
% %     save([save_path '\DATA\' file_data], 'Fs', 'c', 'reso', 'xmin', 'xmax', 'zmin', 'zmax', ...
% %         'CSM', 'freqs', 'Axz', 'scan_plane_Y', 'mic_config', 'flb', 'fub', 'Tblock', 't_start', 't_end', ...
% %         'x', 'z', 'SPLxz', 'dynamic_range', 'maxvalxz');
% % end
% 
% figure;
% contourf(x, y, SPL, (maxval-dynamic_range):1:maxval);
% layoutWingprop(xs, wing, wing_heightTE)
% 
% plot_settings(gca, '$x$ [m]', '$y$ [m]', ['No flow, 1300 [Hz]'], ...
%     [xmin xmax], [ymin ymax], [-0.8 -0.4 0 0.4 0.8], [-0.9 -0.6 -0.3 0 0.3 0.6], ...
%     'on', 'on', 1, [1 (round(maxval)-dynamic_range) round(maxval)], ...
%     '$L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\CB_1300HZ_PROPWINGRPM3']);
% 
% flb = 2635;
% fub = 2685;
% % overlap 50%, blocks of 1 sec start at 0 end at 10
% [CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);
% if remdiag
%     CSM = removeDiagonal(CSM);
% end
% [x, y, A] = FastBeamforming3(CSM, z, freqs, ...
%                 [xmin xmax ymin ymax], reso, mic_config.', c);
% if remdiag
%     A((real(A)<0)) = 1e-20;
% end
% SPL = 20*log10( sqrt(real(A)) / 2e-5 );
% 
% dynamic_range = 12;
% maxval = max(SPL(:));
% figure;
% contourf(x, y, SPL, (maxval-dynamic_range):1:maxval);
% layoutWingprop(xs, wing, wing_heightTE)
% 
% plot_settings(gca, '$x$ [m]', '$y$ [m]', ['No flow, 2660 [Hz]'], ...
%     [xmin xmax], [ymin ymax], [-0.8 -0.4 0 0.4 0.8], [-0.9 -0.6 -0.3 0 0.3 0.6], ...
%     'on', 'on', 1, [1 (round(maxval)-dynamic_range) round(maxval)], ...
%     '$L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\CB_2660HZ_PROPWINGRPM3']);

clear data A SPL;
load(file_path2);
flb = 2635;
fub = 2685;
wing = 0*1;
z = (~wing)*xs(3) + wing*wing_dist;
z_sl = 1.16;
U_c = [0 10 0]*(z-z_sl)/z;

[CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);
if remdiag
    CSM = removeDiagonal(CSM);
end
[x, y, A] = FastBeamforming3(CSM, z, freqs, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
if remdiag
    A((real(A)<0)) = 1e-20;
end
SPL = 20*log10( sqrt(real(A)) / 2e-5 );

dynamic_range = 12;
maxval = max(SPL(:));

figure;
contourf(x, y, SPL, (maxval-dynamic_range):1:maxval);
layoutWingprop(xs, wing, wing_heightTE)

plot_settings(gca, '$x$ [m]', '$y$ [m]', ['No flow, 2660 [Hz]'], ...
    [xmin xmax], [ymin ymax], [-0.8 -0.4 0 0.4 0.8], [-0.9 -0.6 -0.3 0 0.3 0.6], ...
    'on', 'on', 1, [1 (round(maxval)-dynamic_range) round(maxval)], ...
    '$L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\CB_2660HZ_PROPRPM3']);

% 
% source_dist = 1.45;
% source_height = 1.51;
% 
% wing_dist = 1.071;
% wing_heightTE = 1.67;
% 
% xs = [0; source_height-1.57; source_dist];
% clear data A SPL;
% load(file_path3);
% 
% wing = 1;
% z = (~wing)*xs(3) + wing*wing_dist;
% z_sl = 1.16;
% U_c = [0 10 0]*(z-z_sl)/z;
% 
% flb = 1275;
% fub = 1325;
% [CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);
% if remdiag
%     CSM = removeDiagonal(CSM);
% end
% [x, y, A] = FastBeamforming3Conv(CSM, z , freqs, ...
%         [xmin xmax ymin ymax], reso, mic_config.', c, U_c);
% if remdiag
%     A((real(A)<0)) = 1e-20;
% end
% SPL = 20*log10( sqrt(real(A)) / 2e-5 );
% 
% dynamic_range = 12;
% maxval = max(SPL(:));
% 
% figure;
% contourf(x, y, SPL, (maxval-dynamic_range):1:maxval);
% layoutWingprop(xs, wing, wing_heightTE)
% 
% plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Flow, 1300 [Hz]'], ...
%     [xmin xmax], [ymin ymax], [-0.8 -0.4 0 0.4 0.8], [-0.9 -0.6 -0.3 0 0.3 0.6], ...
%     'on', 'on', 1, [1 (round(maxval)-dynamic_range) round(maxval)], ...
%     '$L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\CBFLOW_1300HZ_PROPWINGRPM3']);
% 
% flb = 2635;
% fub = 2685;
% 
% [CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);
% if remdiag
%     CSM = removeDiagonal(CSM);
% end
% [x, y, A] = FastBeamforming3Conv(CSM, z , freqs, ...
%         [xmin xmax ymin ymax], reso, mic_config.', c, U_c);
% if remdiag
%     A((real(A)<0)) = 1e-20;
% end
% SPL = 20*log10( sqrt(real(A)) / 2e-5 );
% 
% dynamic_range = 12;
% maxval = max(SPL(:));
% 
% figure;
% contourf(x, y, SPL, (maxval-dynamic_range):1:maxval);
% layoutWingprop(xs, wing, wing_heightTE)
% 
% plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Flow, 2660 [Hz]'], ...
%     [xmin xmax], [ymin ymax], [-0.8 -0.4 0 0.4 0.8], [-0.9 -0.6 -0.3 0 0.3 0.6], ...
%     'on', 'on', 1, [1 (round(maxval)-dynamic_range) round(maxval)], ...
%     '$L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\CBFLOW_2660HZ_PROPWINGRPM3']);
% 
% clear data A SPL;
% flb = 2635;
% fub = 2685;
% load(file_path4);
% 
% wing = 0;
% z = (~wing)*xs(3) + wing*wing_dist;
% z_sl = 1.16;
% U_c = [0 10 0]*(z-z_sl)/z;
% 
% [CSM, freqs] = developCSM(data, flb, fub, Fs, Tblock, .5, t_start, t_end);
% if remdiag
%     CSM = removeDiagonal(CSM);
% end
% [x, y, A] = FastBeamforming3Conv(CSM, z , freqs, ...
%         [xmin xmax ymin ymax], reso, mic_config.', c, U_c);
% if remdiag
%     A((real(A)<0)) = 1e-20;
% end
% SPL = 20*log10( sqrt(real(A)) / 2e-5 );
% 
% dynamic_range = 12;
% maxval = max(SPL(:));
% 
% figure;
% contourf(x, y, SPL, (maxval-dynamic_range):1:maxval);
% layoutWingprop(xs, wing, wing_heightTE)
% 
% plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Flow, 2660 [Hz]'], ...
%     [xmin xmax], [ymin ymax], [-0.8 -0.4 0 0.4 0.8], [-0.9 -0.6 -0.3 0 0.3 0.6], ...
%     'on', 'on', 1, [1 (round(maxval)-dynamic_range) round(maxval)], ...
%     '$L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\CBFLOW_2660HZ_PROPRPM3']);


function layoutWingprop(xs, wing, wing_heightTE)
xw = [-.615 wing_heightTE-0.245-1.572; .615 wing_heightTE-0.245-1.572; ...
       .615 wing_heightTE-1.572; -.615 wing_heightTE-1.572; ...
      -.615 wing_heightTE-0.245-1.572];
  
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');

patch([-0.147; 0.147; 0.147;-0.147], ...
      [xs(2)-0.005; xs(2)-0.005; xs(2)+0.005; xs(2)+0.005], ...
      0*[1 1 1], 'FaceAlpha', 1);
patch([.3; .3;.359;.359;-.359;-.359;-.3; -.3], ...
      [ylim(1); -.64; -.64; -.635; -.635; -.64; -.64;ylim(1)], ...
      0.25*[1+0.25 1 1], 'FaceAlpha', 0.6);

patch([-.359; -.319; -.319; -.359],[-.64; -.64; .37; .37], 'w', 'FaceAlpha', 0.6);
patch([.319; .359; .359; .319],[-.64; -.64; .37; .37], 'w', 'FaceAlpha', 0.6);
patch([-0.5; 0.5; 0.5; -0.5],[0.2325; 0.2325; 0.2825; 0.2825], 'w', 'FaceAlpha', 0.6);
patch([-0.063/2; 0.063/2; 0.063/2; -0.063/2],[xs(2)+0.021; xs(2)+0.021; xs(2)+.521; xs(2)+.521], 'w', 'FaceAlpha', 0.6);

patch([-0.063/2; 0.063/2; 0.063/4; 0; -0.063/4], ...
      [xs(2)+0.0126; xs(2)+0.0126; xs(2)-0.035; xs(2)+0.0126-0.0616; xs(2)-0.035], ...
      'w', 'FaceAlpha', 0.6);

if wing
    patch(xw(1:end-1,1),xw(1:end-1,2), 0.85*[1 1 1], 'FaceAlpha', 0.75);
end