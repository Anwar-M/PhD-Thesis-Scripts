addpath('.\bf files\');
clear U;
save_data = 0;
save_path = ['.\'];
overlay_bf = 1;

source_dist = 3.00;
source_height = 1.515;
plate_dist = source_dist - 0.351;
plate_height = 1.5;

%%%%%%%%%%%% approximate source position
%xs = [0; -0.057; -2.305+2/100];

xs = [0; source_height-1.572; -source_dist];
%xp = [-15.5 -30; 15.5 -30; 15.5 30; -15.5 30; -15.5 -30]./200;
 xp = [-123 -24.5;123 -24.5;123 24.5;-123 24.5;-123 -24.5]./200;
xp(:,2) = xp(:,2) + plate_height - 1.572;

%%%%%%%%%%%%

mic_config(:,1) = -P_obs(:,1);
mic_config(:,2) = P_obs(:,2);
mic_config(:,3) = P_obs(:,3);

N_mics = size(mic_config, 1);
ps0_freq = ps0_freq';

%CSM = (p_i.' * conj(p_i));
for ff =1:n_freq-1
	clear CSM_aux;
	CSM_aux = (ps0_freq(ff,:).'* conj(ps0_freq(ff,:)));
	for tt=1:64
		for tt2=1:64
			CSM(tt2,tt,ff) = CSM_aux(tt,tt2);
		end
	end
end

%CSM = (ps0_total(1,:).' * conj(ps0_total(1,:)));
freq2=freq(1:end-1);

%% Beamforming
xmin = -1; xmax = 1;
ymin = -1; ymax = 1;
z = 0;

reso = 0.01;
tic;

% Beamforming method, with convection of sound
z_sl = 1.16;
U_c = [0 10 0]*(z-z_sl)/z;

[x, y, A] = FastBeamforming3Conv(CSM, z, freq2,[xmin xmax ymin ymax], reso, mic_config.', c, U_c);

SPL = 20*log10( sqrt(real(A)) / 2e-5 );
disp(toc);

SPL(:) = SPL(:) -  max(SPL(:));
%%

dynamic_range = 12;

maxval = max(SPL(:));
figure;
% imagesc(x, y, SPL, [round(maxval)-dynamic_range round(maxval)]);
contourf(x, y, SPL, (round(maxval)-dynamic_range):1.5:round(maxval));
colormap(jet);

[dummy, vert_ind] = max(SPL);
[~, hor_ind] = max(dummy);
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
placement_x2 = 2*(xlim(2)-xlim(1))/3 + xlim(1);
placement_y2 = 1*(ylim(2)-ylim(1))/30 + ylim(1);
patch([placement_x2-0.025; placement_x2+0.6; placement_x2+0.6; placement_x2-0.025], ...
      [placement_y2-0.1; placement_y2-0.1; placement_y2+0.1; placement_y2+0.1], ...
      [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
text(placement_x2, placement_y2, ['z$_{bf}$ =' num2str(abs(z+plate_dist)) ' [m]'], ...
    'Color','k','FontSize',12,'Interpreter','LaTex');

line([xlim(1); xlim(2)], [xs(2); xs(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([xs(1); xs(1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
patch(xp(1:end-1,1),xp(1:end-1,2), 0.85*[1 1 1], 'FaceAlpha', 0.75);

plot_settings(gca, '$x$ [m]', '$y$ [m]', ['$f$ = ' num2str(freq(1)) ' [Hz]'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
    'on', 'on', 1, [1 round(maxval)-dynamic_range round(maxval)], ...
    ['$\delta$ [dB]'], save_data, [save_path num2str(freq(1)) 'HzPred']);

