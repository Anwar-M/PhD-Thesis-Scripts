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

fileName = 'input_data.txt';

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
    freq = 1950:10:2050;
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
%angle = 0:10:360;
%P_obs(:,1) =1.9.*cos(angle.*2*pi/360);
%P_obs(:,1) = 1.9*(rand(64,1)-0.5);
%P_obs(:,2) = 1.9.*sin(angle.*2*pi/360);
%P_obs(:,2) = 1.9*(rand(64,1)-0.5);
%P_obs(:,3) = -1.1;

%P_obs(2,1) =0;
%P_obs(2,2) = 0;
%P_obs(2,3) = -ac_pos(3);

%P_obs = dlmread('mics_pos.txt');
load('O:\V-Tunnel 20-11 Wing-Source\mic_poses_optim.mat');
P_obs(:,1) = -mic_poses(1,:);
P_obs(:,2) = mic_poses(2,:);
P_obs(:,3) = -1.1;

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
