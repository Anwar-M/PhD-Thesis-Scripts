function EnnesVSPieter

clearvars; close all;

addpath('O:\MATLAB Signal Processing Files');
dirpath = 'O:\PhD Thesis\RESULTS\CH2\ENNESVSPIETER';
save_data = 1;

c = 343.2;
bf_freq = 1000;
x_range = [-1 1];
y_range = [-1 1];
z_bf = 1.0;
y_bf = 0;
z_range = [-1 1] + z_bf;
res = 0.01;
dynamic_range = 6;
stepc = 0.5;

tiz = {'I';'II';'III';'IV'};

load('O:\V-Tunnel 12-12 Wing-Prop\mic_poses_optim.mat');
mic_pos = mic_poses.';
% source_info = [-0.75 0.75 z_bf bf_freq 100];

source_info = [-0.75 0.75 z_bf bf_freq 100; ...
               0 0 z_bf bf_freq 100; ...
               .65 -0.5 z_bf bf_freq 100];

[p, Fs] = simulateArraydata(source_info, mic_pos, c);
[CSM, freqs] = developCSM(p.', bf_freq-5, bf_freq+5, Fs, size(p,2)/Fs, 0);

[X3, Y3, BE] = FastBeamforming3(CSM, z_bf, freqs, [x_range y_range], ...
                             res, mic_pos.', c);
[X3, Y3, BP] = FastBeamforming3mod(CSM, z_bf, freqs, [x_range y_range], ...
                             res, mic_pos.', c);

for I = 1:1
    figure;
	SPLP = 20*log10(sqrt(real(BP))/(4*pi*z_bf*2e-5));
    maxSPLP = 100;
	contourf(X3,Y3,SPLP,[(maxSPLP-dynamic_range):stepc:maxSPLP]);  %imagesc(X3,Z3,B3Lz);
    colormap('hot');
    hold on; 
    scatter(source_info(:,1), source_info(:,2), ...
        'kx', 'LineWidth', 1.5); 
    hold off;
%     maxOut(SPLP,res);
    plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Reference distance'], ...
        x_range, y_range, linspace(x_range(1), x_range(2), 5), ...
        linspace(y_range(1), y_range(2), 5), 16, 'on', 'on', 1, ...
        [1 (maxSPLP-dynamic_range):1:maxSPLP], ...
        '$L_{\mathrm{p}}$ [dB]', save_data, ...
        [dirpath '\PieterForm3']);
    
    figure;
	SPLE = 20*log10(sqrt(real(BE))/2e-5);
    maxSPLE = 100;
	contourf(X3,Y3,SPLE,[(maxSPLE-dynamic_range):stepc:maxSPLE]);  %imagesc(X3,Z3,B3Lz);      
    colormap('hot');
    hold on; 
    scatter(source_info(:,1), source_info(:,2), ...
        'kx', 'LineWidth', 1.5); 
    hold off;
%     maxOut(SPLE,res);
    plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Array center'], ...
        x_range, y_range, linspace(x_range(1), x_range(2), 5), ...
        linspace(y_range(1), y_range(2), 5), 16, 'on', 'on', 1, ...
        [1 (maxSPLE-dynamic_range):1:maxSPLE], ...
        '$L_{\mathrm{p}}$ [dB]', save_data, ...
        [dirpath '\EnnesForm3']);
end

function maxOut(SPL,reso)
[dummy, vert_ind] = max(SPL);
[~, hor_ind] = max(dummy);
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
Xsh = xlim(1) + (hor_ind-1)*reso;
Zsh = ylim(1) + (vert_ind(hor_ind)-1)*reso;

line([xlim(1); xlim(2)], [Zsh; Zsh], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([Xsh; Xsh], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
