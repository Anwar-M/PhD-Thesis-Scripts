function CoherencyBeamforming
clearvars; close all;
save_data = 1;

xs = [-0.12 0.02 1.93; 0.19 0.02 1.93];

pref = 2e-5;
fs = 50e3;

addpath('O:\MATLAB Signal Processing Files');

file_path = ['O:\V-Tunnel Coherent Study\DATA MAT\combi_23-26.mat'];
save_path = ['O:\PhD Thesis\RESULTS\CH2\COHERENCY'];

load(file_path);

load('O:\V-Tunnel Coherent Study\\mic_poses_optim.mat');
mic_config = mic_poses.'; clear mic_poses;
mic_config(:,3) = 0.02;

N_mics = size(mic_config, 1);

%% CSM

c = 345;
t_start = 0;
t_end = 1;
t_block = 1;
t_overlap = 0.5;

f_low = 3000;
f_high = 3001;

% overlap 50%, blocks of 1 sec start at 0 end at 10
[CSMin, freqsin] = developCSM(dataincoh, f_low, f_high, fs, t_block, t_overlap, t_start, t_end);
n_freqsin = length(freqsin);
[CSMc, freqsc] = developCSM(datacoh, f_low, f_high, fs, t_block, t_overlap, t_start, t_end);
n_freqsc = length(freqsc);

%% Beamforming
xmin = -.4; xmax = .4;
ymin = -.4; ymax = .4;
z = 1.93;

reso = 0.01;

tic;
A = [];
for I = 1:length(z)
    [~, ~, Ain] = FastBeamforming3(CSMin, z, freqsin, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
    [x, y, Ac] = FastBeamforming3(CSMc, z, freqsc, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
end
disp(toc);

dynamic_range = 12;

SPLin = 20*log10( sqrt(real(Ain)) / 2e-5 );
SPLc = 20*log10( sqrt(real(Ac)) / 2e-5 );

maxvalin = max(SPLin(:));
minvalin = min(SPLin(:));

figure;
contourf(x, y, SPLin, (round(maxvalin)-dynamic_range):1:round(maxvalin));
colormap('hot');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
line([xlim(1); 0*xlim(2)], [xs(1,2); xs(1,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(1,1); xs(1,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([0*xlim(1); xlim(2)], [xs(2,2); xs(2,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(2,1); xs(2,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Incoherent, $N_{\mathrm{B}} = 1$'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 (round(maxvalin)-dynamic_range) round(maxvalin)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\Inc1Block']);

maxvalc= max(SPLc(:));
% SPLc = SPLc - maxvalc; maxvalc = 0;
minvalc = min(SPLc(:));

figure;
contourf(x, y, SPLc, (round(maxvalc)-dynamic_range):1:round(maxvalc));
colormap('hot');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
line([xlim(1); 0*xlim(2)], [xs(1,2); xs(1,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(1,1); xs(1,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([0*xlim(1); xlim(2)], [xs(2,2); xs(2,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(2,1); xs(2,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Coherent, $N_{\mathrm{B}} = 1$'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 (round(maxvalc)-dynamic_range) round(maxvalc)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\Coh1Block']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_start = 0;
t_end = 1;
t_block = 0.01;
t_overlap = 0.5;

f_low = 3000;
f_high = 3001;

% overlap 50%, blocks of 1 sec start at 0 end at 10
[CSMin, freqsin] = developCSM(dataincoh, f_low, f_high, fs, t_block, t_overlap, t_start, t_end);
n_freqsin = length(freqsin);
[CSMc, freqsc] = developCSM(datacoh, f_low, f_high, fs, t_block, t_overlap, t_start, t_end);
n_freqsc = length(freqsc);


%% Beamforming
xmin = -.4; xmax = .4;
ymin = -.4; ymax = .4;
z = 1.93;

reso = 0.01;

tic;
A = [];
for I = 1:length(z)
    [~, ~, Ain] = FastBeamforming3(CSMin, z, freqsin, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
    [x, y, Ac] = FastBeamforming3(CSMc, z, freqsc, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
end
disp(toc);

%%

dynamic_range = 12;

SPLin = 20*log10( sqrt(real(Ain)) / 2e-5 );
SPLc = 20*log10( sqrt(real(Ac)) / 2e-5 );

maxvalin = max(SPLin(:));
minvalin = min(SPLin(:));

figure;
contourf(x, y, SPLin, (round(maxvalin)-dynamic_range):1:round(maxvalin));
colormap('hot');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
line([xlim(1); 0*xlim(2)], [xs(1,2); xs(1,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(1,1); xs(1,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([0*xlim(1); xlim(2)], [xs(2,2); xs(2,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(2,1); xs(2,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Incoherent, $N_{\mathrm{B}} = 100$'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 (round(maxvalin)-dynamic_range) round(maxvalin)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\Inc100Blocks50overlap']);

maxvalc= max(SPLc(:));
% SPLc = SPLc - maxvalc; maxvalc = 0;
minvalc = min(SPLc(:));

figure;
contourf(x, y, SPLc, (round(maxvalc)-dynamic_range):1:round(maxvalc));
colormap('hot');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
line([xlim(1); 0*xlim(2)], [xs(1,2); xs(1,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(1,1); xs(1,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([0*xlim(1); xlim(2)], [xs(2,2); xs(2,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(2,1); xs(2,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Coherent, $N_{\mathrm{B}} = 100$'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 (round(maxvalc)-dynamic_range) round(maxvalc)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\Coh100Blocks50overlap']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_start = 0;
t_end = 10;
t_block = 0.01;
t_overlap = 0.5;

f_low = 1000;
f_high = 10000;

% overlap 50%, blocks of 1 sec start at 0 end at 10
[CSMin, freqsin] = developCSM(dataincoh, f_low, f_high, fs, t_block, t_overlap, t_start, t_end);
n_freqsin = length(freqsin);
[CSMc, freqsc] = developCSM(datacoh, f_low, f_high, fs, t_block, t_overlap, t_start, t_end);
n_freqsc = length(freqsc);


%% Beamforming
xmin = -.4; xmax = .4;
ymin = -.4; ymax = .4;
z = 1.93;

reso = 0.01;

tic;
A = [];
for I = 1:length(z)
    [~, ~, Ain] = FastBeamforming3(CSMin, z, freqsin, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
    [x, y, Ac] = FastBeamforming3(CSMc, z, freqsc, ...
                [xmin xmax ymin ymax], reso, mic_config.', c);
end
disp(toc);

%%

dynamic_range = 12;

SPLin = 20*log10( sqrt(real(Ain)) / 2e-5 );
SPLc = 20*log10( sqrt(real(Ac)) / 2e-5 );

maxvalin = max(SPLin(:));
minvalin = min(SPLin(:));

figure;
contourf(x, y, SPLin, (round(maxvalin)-dynamic_range):1:round(maxvalin));
colormap('hot');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
line([xlim(1); 0*xlim(2)], [xs(1,2); xs(1,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(1,1); xs(1,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([0*xlim(1); xlim(2)], [xs(2,2); xs(2,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(2,1); xs(2,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Incoherent, $N_{\mathrm{B}} = 1000$'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 (round(maxvalin)-dynamic_range) round(maxvalin)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\IncAllFreq']);

maxvalc= max(SPLc(:));
% SPLc = SPLc - maxvalc; maxvalc = 0;
minvalc = min(SPLc(:));

figure;
contourf(x, y, SPLc, (round(maxvalc)-dynamic_range):1:round(maxvalc));
colormap('hot');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
line([xlim(1); 0*xlim(2)], [xs(1,2); xs(1,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(1,1); xs(1,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([0*xlim(1); xlim(2)], [xs(2,2); xs(2,2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);
line([xs(2,1); xs(2,1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1.5);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Coherent, $N_{\mathrm{B}} = 1000$'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), 16, ...
    'on', 'on', 1, [1 (round(maxvalc)-dynamic_range) round(maxvalc)], '$L_{\mathrm{p}}$ [dB]', save_data, [save_path '\CohAllFreq']);


end
