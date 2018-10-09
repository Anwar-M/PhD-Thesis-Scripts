function multi_beamforming_platesource
% plate or only source
plate = 1;

file_name = 'Mes1.mat';

source_dist = 1.72;
plate_dist = source_dist - 0.63;
plate_height = 1.5;

xs = [0; 1.5-1.572; source_dist];
xp = [-15.5 -30; 15.5 -30; 15.5 30; -15.5 30; -15.5 -30]./200;
xp(:,2) = xp(:,2) + plate_height - 1.572;

addpath('O:\MATLAB Signal Processing Files');

file_path = ['O:\V-Tunnel 30-10 Plate-Source\Data MAT\' file_name];

dynamic_range = 12;

mainresultpath = 'O:\PhD Thesis\RESULTS\CH6\OMNISOURCE';
save_data = 0*1;
Fs = 50e3;

load('O:\V-Tunnel 30-10 Plate-Source\mic_poses_optim.mat');
mic_config = mic_poses.'; clear mic_poses;
mic_config(:,3) = 0.02;
N_mics = size(mic_config, 1);

c = 345;

bfresultpath = [mainresultpath '\' file_name(1:end-4) '_BF'];

fcentre = 2000;
fupper = fcentre + 50;
flower = fcentre - 50;

Nf = length(fcentre);
fprintf('\tUsing %d bands...\n', Nf);

reso = 0.01;

fprintf('\tLoading: %s ...\n', file_name);
load(file_path);
z = (~plate)*xs(3) + plate*plate_dist;

reverseStr = '';
for I = 1:Nf
    msg = sprintf('\tEvaluating frequency %d/%d at %.2f Hz...\n', I, Nf, fcentre(I));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    % overlap 50%, blocks of 1 sec start at 0 end at 10
    [CSM, freqs] = developCSM(data, flower(I), fupper(I), Fs, 0.1, .5, 0, 30);
    xmin = -0.5; xmax = 0.5;
    ymin = -0.5; ymax = 0.5;
    
    [x, y, A] = FastBeamforming3(CSM, z, freqs, ...
        [xmin xmax ymin ymax], reso, mic_config.', c);
    SPL = 20*log10( sqrt(real(A)) / 2e-5 );
    SPL = SPL - max(SPL(:));
    maxval = max(SPL(:));
    valrange = [-12 0];
    figure;
    contourf(x, y, SPL, (round(maxval)-dynamic_range):1:round(maxval));
    layout(z, xs, xp, plate);
    plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Experiment CB'], ...
        [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
        'on', 'on', 1, [1 valrange], '$L_{\mathrm{p, rel}}$ [dB]', save_data, [bfresultpath '_' num2str(round(fcentre(I))) 'HzCB']);
    eval(char(save_data*'close'));
    
    [x, y, A] = adaptive_HR_CleanSC_mod(CSM, z, freqs, ...
            [xmin xmax ymin ymax], reso, mic_config.', c, 4);
    SPL = 20*log10( sqrt(real(A)) / 2e-5 );
    SPL = SPL - max(SPL(:));
    maxval = max(SPL(:));
    SPL = SPL - maxval;
    valrange = [-12 0];
    figure;
    contourf(x, y, SPL, (round(maxval)-dynamic_range):1:round(maxval));
    layout(z, xs, xp, plate);
    plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Experiment HR-CLEAN-SC'], ...
        [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
        'on', 'on', 1, [1 valrange], '$L_{\mathrm{p, rel}}$ [dB]', save_data, [bfresultpath '_' num2str(round(fcentre(I))) 'HzHR']);
    eval(char(save_data*'close'));

end

file_name = 'Mes6.mat';
source_dist = 2.03;
plate_dist = source_dist - 0.2;
plate_height = 1.5;

xs = [0; 1.5-1.572; source_dist];
xp = [-15.5 -30; 15.5 -30; 15.5 30; -15.5 30; -15.5 -30]./200;
xp(:,2) = xp(:,2) + plate_height - 1.572;

file_path = ['O:\V-Tunnel 30-10 Plate-Source\Data MAT\' file_name];
dynamic_range = 12;

bfresultpath = [mainresultpath '\' file_name(1:end-4) '_BF'];

fcentre = 4000;
fupper = fcentre + 50;
flower = fcentre - 50;

Nf = length(fcentre);
fprintf('\tUsing %d bands...\n', Nf);

reso = 0.01;

fprintf('\tLoading: %s ...\n', file_name);
load(file_path);
z = (~plate)*xs(3) + plate*plate_dist;

reverseStr = '';
for I = 1:Nf
    msg = sprintf('\tEvaluating frequency %d/%d at %.2f Hz...\n', I, Nf, fcentre(I));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    % overlap 50%, blocks of 1 sec start at 0 end at 10
    [CSM, freqs] = developCSM(data, flower(I), fupper(I), Fs, 0.1, .5, 0, 30);
    xmin = -0.5; xmax = 0.5;
    ymin = -0.5; ymax = 0.5;
    
    [x, y, A] = FastBeamforming3(CSM, z, freqs, ...
        [xmin xmax ymin ymax], reso, mic_config.', c);
    SPL = 20*log10( sqrt(real(A)) / 2e-5 );
    SPL = SPL - max(SPL(:));
    maxval = max(SPL(:));
    valrange = [-12 0];
    figure;
    contourf(x, y, SPL, (round(maxval)-dynamic_range):1:round(maxval));
    layout(z, xs, xp, plate);
    plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Experiment CB'], ...
        [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
        'on', 'on', 1, [1 valrange], '$L_{\mathrm{p, rel}}$ [dB]', save_data, [bfresultpath '_' num2str(round(fcentre(I))) 'HzCB']);
    eval(char(save_data*'close'));
    
    [x, y, A] = adaptive_HR_CleanSC_mod(CSM, z, freqs, ...
            [xmin xmax ymin ymax], reso, mic_config.', c, 4);
    SPL = 20*log10( sqrt(real(A)) / 2e-5 );
    SPL = SPL - max(SPL(:));
    maxval = max(SPL(:));
    SPL = SPL - maxval;
    valrange = [-12 0];
    figure;
    contourf(x, y, SPL, (round(maxval)-dynamic_range):1:round(maxval));
    layout(z, xs, xp, plate);
    plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Experiment HR-CLEAN-SC'], ...
        [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
        'on', 'on', 1, [1 valrange], '$L_{\mathrm{p, rel}}$ [dB]', save_data, [bfresultpath '_' num2str(round(fcentre(I))) 'HzHR']);
    eval(char(save_data*'close'));

end

fprintf([reverseStr, '\tDone!\n']);

function layout(z, xs, xp, plate)

xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
% placement_x2 = 2*(xlim(2)-xlim(1))/3 + xlim(1);
% placement_y2 = 1*(ylim(2)-ylim(1))/30 + ylim(1);
% patch([placement_x2-0.025; placement_x2+0.6; placement_x2+0.6; placement_x2-0.025], ...
%       [placement_y2-0.1; placement_y2-0.1; placement_y2+0.1; placement_y2+0.1], ...
%       [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
% text(placement_x2, placement_y2, ['z$_{bf}$ = ' num2str(abs(z)) ' [m]'], ...
%     'Color','k','FontSize',14,'Interpreter','LaTex');

line([xlim(1); xlim(2)], [xs(2); xs(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([xs(1); xs(1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);

if plate
    % line(xp(:,1),xp(:,2), 'Color', 0.85*[1 1 1], 'LineStyle', '-','LineWidth', 1.5);
    patch(xp(1:end-1,1),xp(1:end-1,2), 0.85*[1 1 1], 'FaceAlpha', 0.75);
end
