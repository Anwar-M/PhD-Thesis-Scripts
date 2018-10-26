function ConvectionAndDiagRemoval

clearvars; close all;

addpath('O:\MATLAB Signal Processing Files');
dirpath = 'O:\PhD Thesis\RESULTS\CH2\RANDOMCSM';
save_data = 0;

c = 343.2;
bf_freq = 2000;
N_grid1D = 100;
x_range = 1*[-1 1];
y_range = 1*[-1 1];
z_bf = 1;

res = 0.01;
dynamic_range = 24;
stepc = 1.5;

tiz = {'Original CSM'; 'Random powers'; 'Random phase I';'Random phase II'};
filzz = {'OrigCSM'; 'RandPowCSM'; 'RandPh1CSM'; 'RandPh2CSM'};

load('O:\V-Tunnel 12-12 Wing-Prop\mic_poses_optim.mat');
mic_pos = mic_poses.';
source_info = [0 0 z_bf bf_freq 100];
% source_info = [-0.5 .6 z_bf bf_freq 100; -0.1 -0.35 z_bf bf_freq 100; 0.6 0.1 z_bf bf_freq 100];
[p, Fs] = simulateArraydata(source_info, mic_pos, c, 50e3, 5);
p = p + wgn(size(p,1),size(p,2), 40);

[CSM, freqs] = developCSM(p.', bf_freq-2, bf_freq+2, Fs, 1, 0.5);
[X, Y, B(:,:,1)] = FastBeamforming3(CSM, z_bf, freqs, [x_range y_range], ...
                             0.01, mic_pos.', c);

CSM = removeDiagonal(CSM);
                         
[X, Y, Bee] = FastBeamforming3(CSM, z_bf, freqs, [x_range y_range], ...
                             0.01, mic_pos.', c);
Bee((real(Bee)<0)) = 1e-20;
B(:,:,2)  = Bee;
% randAmp = (max(abs(CSM(:)))-min(abs(CSM(:))))*rand(64,64)+min(abs(CSM(:)));
% randAmp = triu(randAmp) + triu(randAmp,1).';
% CSM2 = randAmp.*exp(1i*angle(CSM));
% 
% [X, Y, B(:,:,2)] = FastBeamforming3(CSM2, z_bf, freqs, [x_range y_range], ...
%                              0.01, mic_pos.', c);
%                          
% randAng = 2*pi*rand(64,64) - pi;
% randAng = triu(randAng,1) - triu(randAng,1).' ;
% CSM3 = abs(CSM).*exp(1i*randAng);
% 
% [X, Y, B(:,:,3)] = FastBeamforming3(CSM3, z_bf, freqs, [x_range y_range], ...
%                              0.01, mic_pos.', c);
%                          
% randAngLittle = (2*pi*rand(64,64) - pi)*3/4;
% randAngLittle = triu(randAngLittle,1) - triu(randAngLittle,1).' ;
% CSM4 = abs(CSM).*exp(1i*(angle(CSM)+randAngLittle));
% 
% [X, Y, B(:,:,4)] = FastBeamforming3(CSM4, z_bf, freqs, [x_range y_range], ...
%                              0.01, mic_pos.', c);


for I = 1:2
    figure;
	SPL = 20*log10(sqrt(real(abs(B(:,:,I))))/2e-5);
    maxSPL = round(max(SPL(:))*10)/10;
	contourf(X,Y,SPL,[(maxSPL-dynamic_range):stepc:maxSPL]);  %imagesc(X,y,B3Lz);
%     imagesc(X,Y,SPL);
    colormap('hot');
    axis equal; axis([x_range y_range]);
    plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [tiz{I}], ...
        x_range, y_range, linspace(x_range(1), x_range(2), 5), ...
        linspace(y_range(1), y_range(2), 5), 16, 'on', 'on', 1, ...
        [1 (maxSPL-dynamic_range):dynamic_range/6:maxSPL], ...
        '$L_{\mathrm{p}}$ [dB]', save_data, ...
        [dirpath '\1SR' filzz{I}]);
    clear SPL
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
