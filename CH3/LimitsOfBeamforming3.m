function LimitsOfBeamforming3
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
bf_freq = 2000;
N_grid1D = 100;
x_range = 1*[-1 1];
y_range = 1*[-1 1];
z_range = 1;
dBrange = 50;
res = 0.01;
stepc = 1.5;

R_ap = 1;

source_info = [0 0 z_range bf_freq 100];

N_mic = 2000;
[mic_pos, N_mic] = generateArray('reg', N_mic, R_ap, 500);
[p, Fs] = simulateArraydata(source_info, mic_pos, c);
[CSM, freqs] = developCSM(p.', bf_freq-5, bf_freq+5, Fs, size(p,2)/Fs, 0);

eval(evalBeamformSteering(1));
BB1 = 20*log10(sqrt(real(B1(101,:)))/2e-5); clear B1;
eval(evalBeamformSteering(2));
BB2 = 20*log10(sqrt(real(B2(101,:)))/2e-5); clear B2;
eval(evalBeamformSteering(3));
BB3 = 20*log10(sqrt(real(B3(101,:)))/2e-5); clear B3;
eval(evalBeamformSteering(4));
BB4 = 20*log10(sqrt(real(B4(101,:)))/2e-5); clear B4;

figure;
plot(X,BB1,X,BB2,X,BB3,X,BB4);
hl = legend('Formulation I', 'Formulation II','Formulation III','Formulation IV', 'Location', 'NorthEast');
set(hl, 'Interpreter', 'LaTex');
line([-1 1],[82.43 82.43], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
plot_settings_font(gca, '$x$ [m]', '$L_{\mathrm{p}}$ [dB]', ['$y = 0$ [m]'], ...
    x_range, [70 100], linspace(x_range(1), x_range(2), 5), ...
    [70 82.43 90 100], 16, 'on', 'on', 0, ...
    [0], [], save_data, [dirpath '\miccfg_steerings']);
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