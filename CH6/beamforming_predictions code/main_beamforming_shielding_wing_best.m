%----------- Noise Shielding for real AC V02 ----------
%----------- Coded by Ana Vieira 05-09-2016 -----------

%-------------------- Requires ------------------------
%- abs_atm_function.m
%- asym_diff.m 
%- atm_effects.m
%- factor_eta.m
%- find_shadow.m 
%- read_input.m
%-------------------------------------------------------

clear all;
tic
%--- Read input file ---

Z_SOURCE_ARRAY = 2.305;
source_height = 1.515;
freq = 1950:10:2050;

addpath('O:\MATLAB Signal Processing Files');
save_data = 1;
save_path = ['O:\PhD Thesis\RESULTS\CH6\OMNISOURCE\'];

fileName = 'input_data_wing.txt';

[numerical_inp,array_inp,array_inp2,psd_file,geo_file,sources]=read_input(fileName);

p_ref   	= numerical_inp(1);
c 	    	= numerical_inp(2);
ac_pos  	= array_inp(1,:);
hum	    	= numerical_inp(3);
T	    	= numerical_inp(4);
U	    	= array_inp(2,:);
total_t 	= numerical_inp(5);
delta_t 	= numerical_inp(6);
source_N	= sources;
Nx			= numerical_inp(7);
Ny			= numerical_inp(8);
P1			=array_inp2(1,:);
P2			=array_inp2(2,:);
P3			=array_inp2(3,:);
P4			=array_inp2(4,:);
n_sources	=numerical_inp(9);
%----


mach = U/c;					%Mach number

time = 0:delta_t:total_t;	%Time discretization
n_time = size(time,2);		%Number of time steps


for ss=1:n_sources


	%Read source data (PSD vs frequency)
	aux = dlmread(psd_file{ss});
% 	freq = aux(:,1);
% 	PSD = aux(:,2);
    PSD = rand(numel(freq),1);

	n_freq = length(freq);	%Number of frequencies

	%Read shielding body geometry
	pt_aux(ss,:,:)=dlmread(geo_file{ss});	
	
	source_pos(1,ss,:) = source_N(ss,:);
	source_pos(1,ss,1) = source_pos(1,ss,1);% +ac_pos(1)-max(pt_aux(ss,:,1));	%Place source in the right initial position

	pt(ss,1,:,1) = pt_aux(ss,:,1)+ac_pos(1)-max(pt_aux(ss,:,1)); %translate aircraft
	pt(ss,1,:,2) = pt_aux(ss,:,2)+ac_pos(2);
	pt(ss,1,:,3) = pt_aux(ss,:,3);
	
	

	
	%Update of the source position
	for tt=2:n_time
		source_pos(tt,ss,1) = source_pos(tt-1,ss,1);%+U(1)*delta_t;
		source_pos(tt,ss,2) = source_pos(tt-1,ss,2);
		source_pos(tt,ss,3) = source_pos(tt-1,ss,3);
	end


	N_points = size(pt,3);	%Number of points in the shielding body

	%AC position relative to the source
	pt_1(ss,1,:,1) = pt(ss,1,:,1)- source_pos(1,ss,1);
	pt_1(ss,1,:,2) = pt(ss,1,:,2)- source_pos(1,ss,2);
	pt_1(ss,1,:,3) = pt(ss,1,:,3)- source_pos(1,ss,3);


	%Update aircraft position
	for tt=2:n_time
		pt(ss,tt,:,1) = pt(ss,tt-1,:,1);% + U(1)*delta_t;
		pt(ss,tt,:,2) = pt(ss,tt-1,:,2);
		pt(ss,tt,:,3) = pt(ss,tt-1,:,3);
		
		pt_1(ss,tt,:,1) = pt(ss,tt,:,1) - source_pos(tt,ss,1);
		pt_1(ss,tt,:,2) = pt(ss,tt,:,2) - source_pos(tt,ss,2);
		pt_1(ss,tt,:,3) = pt(ss,tt,:,3) - source_pos(tt,ss,3);
	
	end
	

	
end
%--- Define grid of observers ---

P1(3)= -ac_pos(3);
P2(3)= -ac_pos(3);
P3(3)= -ac_pos(3);
P4(3)= -ac_pos(3);

N_obs = Nx*Ny;

dx = (P3(1) - P1(1) )/(Nx-1);
dy = (P2(2) - P1(2) )/(Ny-1);

[X,Y] = meshgrid(P1(1):dx:P4(1),P1(2):dy:P2(2));

n=1;
for k=1:Ny
	for l=1:Nx
		P_obs(n,1:3) = [P1(1)+dx*(l-1) P1(2)+dy*(k-1) -ac_pos(3)];
		n = n+1;
	end
end
clear P_obs;

%P_obs = dlmread('mics_pos.txt');
load('O:\V-Tunnel 20-11 Wing-Source\mic_poses_optim.mat');
P_obs(:,1) = -mic_poses(1,:);
P_obs(:,2) = mic_poses(2,:);
P_obs(:,3) = source_N(1,3)-Z_SOURCE_ARRAY+(2/100);

N_obs = size(P_obs,1);

for ss=1:n_sources
	for tt=1:n_time
		for nn=1:N_obs
			P_obs_1(ss,tt,nn,1) = P_obs(nn,1) - source_pos(tt,ss,1);
			P_obs_1(ss,tt,nn,2) = P_obs(nn,2) - source_pos(tt,ss,2);
			P_obs_1(ss,tt,nn,3) = P_obs(nn,3) - source_pos(tt,ss,3);
		end
	end
end
%----------------------------------

for tt=1:n_time
	
	for ss=1:n_sources
	
		%---- Find shadow zone ----
		pt2 = squeeze(pt_1(ss,:,:,:));
		P_obs2 = squeeze(P_obs_1(ss,:,:,:));
		source_pos_aux(tt,:) = squeeze(source_pos(tt,ss,:));
		
		in = find_shadow(N_points,N_obs,tt,pt2,P_obs2);
		
		%--------------------------

		for ff=1:n_freq-1
			for n=1:N_obs
				
				
				%Include atmospheric effects in the incident pressure (doppler effect and absorption)
				[SPL_obs(n,ff),freq_doppler,k_number,t_obs(n,tt)]= atm_effects(tt,ff,P_obs(n,:),source_pos_aux(tt,:),U,delta_t,c,freq,PSD(ff),T,hum);
	
				%Calculate diffraction
				[I_s,I_end] = asymp_diff(tt,pt2,N_points,P_obs2,n,k_number);
				
				%Distance between source and observer
				no_P_obs = sqrt(P_obs2(tt,n,1)^2 + P_obs2(tt,n,2)^2 + P_obs2(tt,n,3)^2);
				
				%Incident acoustic pressure
				p_i(n) = exp(1i*k_number*no_P_obs)/no_P_obs;
				
				%Apply Babinet's principle
				integral = sum(I_end + I_s);	
				if( in(n)==0)
					integral2 = integral*exp(1i*k_number*no_P_obs);
				else
					integral2 = integral*exp(1i*k_number*no_P_obs)+p_i(n);
				end
				ps0(n)=(p_i(n)-integral2);
				ps0_freq(n,ff) = ps0(n);
				%Calculate Attenuation [dB]
				spl(tt,ss,n,ff) = 20*log10(abs((ps0(n))/p_i(n)));
				
				%if(spl(tt,ss,n,ff)>0)
				%	spl(tt,ss,n,ff) = 10000;
				%end
				
				%fprintf('%d, spl %f\n',n,spl(tt,n,ff));
				
				%SPL [dB] at the observers with attenuation
				%if(spl(tt,ss,n,ff)<0)
					noise_final_1(tt,ss,n,ff) = SPL_obs(n,ff) + spl(tt,ss,n,ff);
				%else
				%	noise_final_1(tt,ss,n,ff) = SPL_obs(n,ff);
				%end
				noise_no_shield(tt,ss,n,ff) = SPL_obs(n,ff);
				
				%SPL [dBA] at the observers with attenuation
				delta_LA = -145.528 + 98.262*log10(freq_doppler(ff)) -19.509*(log10(freq_doppler(ff)))^2 + 0.975*(log10(freq_doppler(ff)))^3;
				
				noise_no_shield_dBA(tt,ss,n,ff) = SPL_obs(n,ff)+delta_LA;
				
				SPL_dBA_1(tt,ss,n,ff) = noise_final_1(tt,ss,n,ff) + delta_LA;
				
				
				%if(spl(tt,ss,n,ff)<0)
					att_dBA(tt,ss,n,ff) = spl(tt,ss,n,ff) + delta_LA;
				%else
				%	att_dBA(tt,ss,n,ff) =delta_LA;
				%end
			end
			fprintf('Freq doppler %f\n',freq_doppler(ff));
		end
		fprintf('\n===================\n');
		fprintf('Time %f\n',tt*delta_t);
		fprintf('===================\n');
		
for tt=1:n_time
	for nn=1:N_obs
		for ff=1:n_freq-1
			ps0_total(nn) = sum(ps0_freq(nn,:));
		end
	end
end

		
		
	end

end

%%

clear U;
overlay_bf = 1;

source_dist = Z_SOURCE_ARRAY;
plate_dist = source_dist - 0.875;
wing_TE = 1.605;

%%%%%%%%%%%% approximate source position
%xs = [0; -0.057; -2.305+2/100];

xs = [0; source_height-1.572; -source_dist];
%xp = [-15.5 -30; 15.5 -30; 15.5 30; -15.5 30; -15.5 -30]./200;
  xp = [-123 -24.5;123 -24.5;123 24.5;-123 24.5;-123 -24.5]./200;
xp(:,2) = xp(:,2) + wing_TE - 1.572 - 0.1225;

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
xmin = -0.7; xmax = 0.7;
ymin = -0.7; ymax = 0.7;
z = 0;

reso = 0.01;
tic;

% Beamforming method, uncomment the one you want. Comment others! nu is needed for Functional Beamforming

[x, y, A] = FastBeamforming3(CSM, z, freq2,[xmin xmax ymin ymax], reso, mic_config.', c);

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
% placement_x2 = 2*(xlim(2)-xlim(1))/3 + xlim(1);
% placement_y2 = 1*(ylim(2)-ylim(1))/30 + ylim(1);
% patch([placement_x2-0.025; placement_x2+0.6; placement_x2+0.6; placement_x2-0.025], ...
%       [placement_y2-0.1; placement_y2-0.1; placement_y2+0.1; placement_y2+0.1], ...
%       [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
% text(placement_x2, placement_y2, ['z$_{bf}$ =' num2str(abs(z+plate_dist)) ' [m]'], ...
%     'Color','k','FontSize',12,'Interpreter','LaTex');

line([xlim(1); xlim(2)], [xs(2); xs(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([xs(1); xs(1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
patch(xp(1:end-1,1),xp(1:end-1,2), 0.85*[1 1 1], 'FaceAlpha', 0.75);

plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Model CB'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
    'on', 'on', 1, [1 round(maxval)-dynamic_range round(maxval)], ...
    ['$L_{\mathrm{p,rel}}$ [dB]'], save_data, [save_path 'PR3WING_BF_2000HzCB']);

%%

[x, y, A] = adaptive_HR_CleanSC_mod(CSM, z, freq(1:end-1),[xmin xmax ymin ymax], reso, mic_config.', c, 6);
SPL = 20*log10( sqrt(real(A)) / 2e-5 );
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
% placement_x2 = 2*(xlim(2)-xlim(1))/3 + xlim(1);
% placement_y2 = 1*(ylim(2)-ylim(1))/30 + ylim(1);
% patch([placement_x2-0.025; placement_x2+0.6; placement_x2+0.6; placement_x2-0.025], ...
%       [placement_y2-0.1; placement_y2-0.1; placement_y2+0.1; placement_y2+0.1], ...
%       [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
% text(placement_x2, placement_y2, ['z$_{bf}$ =' num2str(abs(z+plate_dist)) ' [m]'], ...
%     'Color','k','FontSize',12,'Interpreter','LaTex');

line([xlim(1); xlim(2)], [xs(2); xs(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
line([xs(1); xs(1)], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
patch(xp(1:end-1,1),xp(1:end-1,2), 0.85*[1 1 1], 'FaceAlpha', 0.75);

plot_settings(gca, '$x$ [m]', '$y$ [m]', ['Model HR-CLEAN-SC'], ...
    [xmin xmax], [ymin ymax], linspace(xmin, xmax, 5), linspace(ymin, ymax, 5), ...
    'on', 'on', 1, [1 round(maxval)-dynamic_range round(maxval)], ...
    ['$L_{\mathrm{p,rel}}$ [dB]'], save_data, [save_path 'PR3WING_BF_2000HzHR']);

