function SteeringComparisons

clearvars; close all;

addpath('O:\MATLAB Signal Processing Files');
dirpath = 'O:\PhD Thesis\RESULTS\CH2\STEERING';
save_data = 0;

c = 343.2;
bf_freq = 1200;
x_range = 1/2*[-1 1];
y_range = 1/2*[-1 1];
z_bf = 1.0;
y_bf = 0;
z_range = 1/2*[-1 1] + z_bf;
res = 0.01;
dynamic_range = 50;
stepc = 1.5;

tiz = {'I';'II';'III';'IV'};

load('O:\V-Tunnel 12-12 Wing-Prop\mic_poses_optim.mat');
mic_pos = mic_poses.';
source_info = [-0.25 0.25 z_bf bf_freq 100];

% source_info = [-0.35 0 z_bf bf_freq 100; ...
%                -0.05 0 z_bf bf_freq 100; ...
%                0.25 0 z_bf bf_freq 100];

[p, Fs] = simulateArraydata(source_info, mic_pos, c);
[CSM, freqs] = developCSM(p.', bf_freq-5, bf_freq+5, Fs, size(p,2)/Fs, 0);

[X1, Y1, B(:,:,1)] = FastBeamforming1(CSM, z_bf, freqs, [x_range y_range], ...
                             res, mic_pos.', c);
[X2, Y2, B(:,:,2)] = FastBeamforming2(CSM, z_bf, freqs, [x_range y_range], ...
                             res, mic_pos.', c);
[X3, Y3, B(:,:,3)] = FastBeamforming3(CSM, z_bf, freqs, [x_range y_range], ...
                             res, mic_pos.', c);
[X4, Y4, B(:,:,4)] = FastBeamforming4(CSM, z_bf, freqs, [x_range y_range], ...
                             res, mic_pos.', c);
                         
[X1, Z1, Bz(:,:,1)] = FastBeamforming1(CSM, y_bf, freqs, [x_range z_range], ...
                             res, [mic_pos(:,1) mic_pos(:,3) mic_pos(:,2)].', c);
[X2, Z2, Bz(:,:,2)] = FastBeamforming2(CSM, y_bf, freqs, [x_range z_range], ...
                             res, [mic_pos(:,1) mic_pos(:,3) mic_pos(:,2)].', c);
[X3, Z3, Bz(:,:,3)] = FastBeamforming3(CSM, y_bf, freqs, [x_range z_range], ...
                             res, [mic_pos(:,1) mic_pos(:,3) mic_pos(:,2)].', c);
[X4, Z4, Bz(:,:,4)] = FastBeamforming4(CSM, y_bf, freqs, [x_range z_range], ...
                             res, [mic_pos(:,1) mic_pos(:,3) mic_pos(:,2)].', c);

for I = 1:4
    figure;
	SPL = 20*log10(sqrt(real(Bz(:,:,I)))/2e-5);
    maxSPL = round(max(SPL(:))*10)/10;
	contourf(X3,Z3,SPL,[(maxSPL-dynamic_range):stepc:maxSPL]);  %imagesc(X3,Z3,B3Lz);
    colormap('hot');
    hold on; 
    scatter(source_info(:,1), source_info(:,3), ...
        'kx', 'LineWidth', 1.5); 
    hold off;
    maxOut(SPL,res);
    axis equal; axis([x_range z_range]);
    plot_settings_font(gca, '$x$ [m]', '$z$ [m]', ['Formulation ' tiz{I}], ...
        x_range, z_range, linspace(x_range(1), x_range(2), 5), ...
        linspace(z_range(1), z_range(2), 5), 16, 'on', 'on', 1, ...
        [1 (maxSPL-dynamic_range):10:maxSPL], ...
        '$L_{\mathrm{p}}$ [dB]', save_data, ...
        [dirpath '\Form' tiz{I} 'xz']);
    
    figure;
	SPL = 20*log10(sqrt(real(B(:,:,I)))/2e-5);
    maxSPL = round(max(SPL(:))*10)/10;
	contourf(X3,Y3,SPL,[(maxSPL-dynamic_range):stepc:maxSPL]);  %imagesc(X3,Z3,B3Lz);      
    colormap('hot');
    hold on; 
    scatter(source_info(:,1), source_info(:,2), ...
        'kx', 'LineWidth', 1.5); 
    hold off;
    maxOut(SPL,res);
    axis equal; axis([x_range y_range]);
    plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Formulation ' tiz{I}], ...
        x_range, y_range, linspace(x_range(1), x_range(2), 5), ...
        linspace(y_range(1), y_range(2), 5), 16, 'on', 'on', 1, ...
        [1 (maxSPL-dynamic_range):10:maxSPL], ...
        '$L_{\mathrm{p}}$ [dB]', save_data, ...
        [dirpath '\Form' tiz{I} 'xy']);
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
