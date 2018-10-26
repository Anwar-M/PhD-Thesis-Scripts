function LimitsOfBeamforming2
% Testing the limits of beamforming and mic arrays
%

%   Anwar Malgoezar, October 2018. 
%   Group ANCE

clearvars; close all;

addpath('O:\MATLAB Signal Processing Files');
dirpath = 'O:\PhD Thesis\RESULTS\CH3\App';
save_data = 5;
size_le_marker = 10;
option = 'reg';

c = 343.2;
bf_freq = 1000;
N_grid1D = 100;
x_range = 1*[-1 1];
y_range = [-1 1];
z_range = 1;
dBrange = 50;
res = 0.01;
stepc = 1.5;

N_mic = 500;
R_ap1 = 0.5;
[mic_pos1, N_mic] = generateArray('circ', N_mic, R_ap1, 75);
R_ap2 = 1;
[mic_pos2, N_mic] = generateArray('circ', N_mic, R_ap2, 75);

figure; plot(mic_pos1(:,1),mic_pos1(:,2),'b.',mic_pos2(:,1),mic_pos2(:,2),'r.','MarkerSize',size_le_marker);
hl = legend(['Configuration I'], ['Configuration II'], 'Location', 'NorthEast');
set(hl, 'Interpreter', 'LaTex');
plot_settings_font(gca, '$x$ [m]', '$y$ [m]', ['Microphones'], ...
    1.1*x_range, 1.1*y_range, linspace(x_range(1), x_range(2), 5), ...
    linspace(y_range(1), y_range(2), 5), 16, 'on', 'on', 1, ...
    [0], [], save_data, [dirpath '\miccfg_twocircs']);
close;

source_info = [0 0 z_range bf_freq 100];
[p, Fs] = simulateArraydata(source_info, mic_pos1, c);
[CSM, freqs] = developCSM(p.', bf_freq-5, bf_freq+5, Fs, size(p,2)/Fs, 0);

[X, Y, B3] = FastBeamforming3(CSM, z_range, freqs, [x_range y_range], res, mic_pos1.', c);

BB1 = 20*log10(sqrt(real(B3(101,:)))/2e-5); clear B3;

source_info = [0 0 z_range bf_freq 100];
[p, Fs] = simulateArraydata(source_info, mic_pos2, c);
[CSM, freqs] = developCSM(p.', bf_freq-5, bf_freq+5, Fs, size(p,2)/Fs, 0);

[X, Y, B3] = FastBeamforming3(CSM, z_range, freqs, [x_range y_range], res, mic_pos2.', c);

BB2 = 20*log10(sqrt(real(B3(101,:)))/2e-5); clear B3;



figure;
plot(X,BB1,X,BB2);
hl = legend('Small aperture', 'Large aperture', 'Location', 'NorthEast');
set(hl, 'Interpreter', 'LaTex');
% line([-1 1],[82.43 82.43], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
plot_settings_font(gca, '$x$ [m]', '$L_{\mathrm{p}}$ [dB]', ['$y = 0$ [m]'], ...
    x_range, [70 100], linspace(x_range(1), x_range(2), 5), ...
    [70 80 90 100], 16, 'on', 'on', 0, ...
    [0], [], save_data, [dirpath '\miccfg_y0twocircs']);
close;

function [mic_pos, N_mic] = generateArray(theway, N_mic, R_ap, N_circ)

switch theway
    case 'reg'
        micperdim = ceil(sqrt(N_mic));
        mic_pos = ...
            [repmat(linspace(-R_ap, R_ap, micperdim).',micperdim,1) ...
            reshape(ones(micperdim,1)*linspace(-R_ap, R_ap, micperdim), micperdim^2, 1)];

        mic_pos(sqrt(mic_pos(:,1).^2+mic_pos(:,2).^2)>R_ap, :) = [];

        [xe, ye] = pol2cart(linspace(0,2*pi,N_circ).',R_ap*ones(N_circ,1));
        mic_pos(end:end+N_circ-1, 1:2) = [xe ye];
        N_mic = size(mic_pos,1);
    case 'rand'
        mic_pos(:,1:2) = R_ap*(2*rand(N_mic,2)-1);
        mic_pos(sqrt(mic_pos(:,1).^2+mic_pos(:,2).^2)>R_ap, :) = [];

        [xe, ye] = pol2cart(linspace(0,2*pi,100).',R_ap*ones(100,1));
        mic_pos(end:end+99, 1:2) = [xe ye];
        N_mic = size(mic_pos,1);
    case 'circ'
        TH = 2*pi*rand(1,N_mic);
        RA = R_ap*ones(1,N_mic);
        [mic_pos(:,1), mic_pos(:,2)] = pol2cart(TH,RA);
end
mic_pos(:,3) = 0;

function str = evalBeamformSteering(num)
switch num
    case 1
        str ='[X, Y, B1] = FastBeamforming1(CSM, z_range, freqs, [x_range y_range],res, mic_pos.'', c);';
    case 2
        str = '[X, Y, B2] = FastBeamforming2(CSM, z_range, freqs, [x_range y_range],res, mic_pos.'', c);';
    case 3
        str = '[X, Y, B3] = FastBeamforming3(CSM, z_range, freqs, [x_range y_range],res, mic_pos.'', c);';
    case 4
        str = '[X, Y, B4] = FastBeamforming4(CSM, z_range, freqs, [x_range y_range],res, mic_pos.'', c);';
    otherwise
        error('One to four!');
end