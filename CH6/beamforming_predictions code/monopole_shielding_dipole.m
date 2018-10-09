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
	freq = aux(:,1);
	PSD = aux(:,2);

	n_freq = size(freq,1);	%Number of frequencies

	%Read shielding body geometry
	pt_aux(ss,:,:)=dlmread(geo_file{ss});	
	
	source_pos(1,ss,:) = source_N(ss,:);
	source_pos(1,ss,1) = source_pos(1,ss,1);% +ac_pos(1)-max(pt_aux(ss,:,1));	%Place source in the right initial position

	pt(ss,1,:,1) = pt_aux(ss,:,1)+ac_pos(1)-max(pt_aux(ss,:,1)); %translate aircraft
	pt(ss,1,:,3) = pt_aux(ss,:,2)+ac_pos(2);
	pt(ss,1,:,2) = pt_aux(ss,:,3);
	
	

	
	%Update of the source position
	for tt=2:n_time
		source_pos(tt,ss,1) = source_pos(tt-1,ss,1);%+U(1)*delta_t;
		source_pos(tt,ss,3) = source_pos(tt-1,ss,3);
		source_pos(tt,ss,2) = source_pos(tt-1,ss,2);
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

%[X,Y] = meshgrid(P1(1):dx:P4(1),P1(2):dy:P2(2));

%n=1;
%for k=1:Ny
%	for l=1:Nx
%		P_obs(n,1:3) = [P1(1)+dx*(l-1) P1(2)+dy*(k-1) -ac_pos(3)];
%		n = n+1;
%	end
%end
clear P_obs;
%[theta,phi] = ndgrid(0 : pi/20 :2*pi , 0 : pi/20 : 2*pi);
%R = 0.1;

% Z = R.*cos(phi).*sin(theta);
% Y = R.*sin(phi).*sin(theta) ;
% X = R.*cos(theta);
%Z = R.*cos(theta);
%Y = R.*sin(phi).*sin(theta) ;
%X = R.*cos(phi).*sin(theta);

%aux =1;
%for nn=1:size(X,1)%%%
%	for mm=1:size(X,2)
%		P_obs(aux,1) = X(nn,mm);
%		P_obs(aux,2) = Y(nn,mm);
%		P_obs(aux,3) = Z(nn,mm);
%		aux = aux+1;
%	end
%end

%N_obs = size(P_obs,1);
%[X,Y,Z] = meshgrid(P_obs(:,1),P_obs(:,2),P_obs(:,3));

%P_obs(:,1) =-2:0.01:2;
%P_obs(:,2) = 0;
%P_obs(:,3) = -ac_pos(3);

%P_obs(2,1) =0;
%P_obs(2,2) = 0;
%P_obs(2,3) = -ac_pos(3);

%P_obs = dlmread('mics_pos.txt');
load('O:\Wing Meas - 20_21Nov\mic_poses_optim.mat');
P_obs(:,1) = -mic_poses(1,:);
P_obs(:,3) = mic_poses(2,:);
P_obs(:,2) = -0.92;

N_obs = size(P_obs,1);

%P_obs_aux = P_obs*rotz(90);
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
				
				%Distance between source and observer
				no_P_obs = sqrt(P_obs2(tt,n,1)^2 + P_obs2(tt,n,2)^2 + P_obs2(tt,n,3)^2);
				
				%Include atmospheric effects in the incident pressure (doppler effect and absorption)
				[SPL_obs(n,ff),freq_doppler,k_number,t_obs(n,tt)]= atm_effects(tt,ff,P_obs(n,:),source_pos_aux(tt,:),U,delta_t,c,freq,PSD(ff),T,hum);
	
				%Calculate diffraction
				[I_s,I_end] = asymp_diff_dipole(tt,pt2,N_points,P_obs2,n,k_number,no_P_obs);
				
				%Distance between source and observer
				no_P_obs = sqrt(P_obs2(tt,n,1)^2 + P_obs2(tt,n,2)^2 + P_obs2(tt,n,3)^2);
	
				%Incident acoustic pressure
				p_i1(n) = (-dot([0,0,1],squeeze(P_obs2(tt,n,:)))*(1-1i*k_number*no_P_obs)/no_P_obs^2)*exp(1i*k_number*no_P_obs)/no_P_obs;
				
				
				p_i2(n) = (-dot([0.8,0,0.6],squeeze(P_obs2(tt,n,:)))*(1-1i*k_number*no_P_obs)/no_P_obs^2)*exp(1i*k_number*no_P_obs)/no_P_obs;
				
				p_i3(n) = (-dot([0.6,0,0.8],squeeze(P_obs2(tt,n,:)))*(1-1i*k_number*no_P_obs)/no_P_obs^2)*exp(1i*k_number*no_P_obs)/no_P_obs;
				
				p_i4(n) = (-dot([0.4,0,0.91652],squeeze(P_obs2(tt,n,:)))*(1-1i*k_number*no_P_obs)/no_P_obs^2)*exp(1i*k_number*no_P_obs)/no_P_obs;
				
				p_i5(n) = (-dot([0.2,0,0.979796],squeeze(P_obs2(tt,n,:)))*(1-1i*k_number*no_P_obs)/no_P_obs^2)*exp(1i*k_number*no_P_obs)/no_P_obs;
				
				p_i6(n) = (dot([0,-0.707106,0.707106],squeeze(P_obs2(tt,n,:)))*(1-1i*k_number*no_P_obs)/no_P_obs^2)*exp(1i*k_number*no_P_obs)/no_P_obs;
				
				%p_i7(n) =  k_number^2*1*sin(phi(n)).*exp(-1i*k_number*no_P_obs)./(4*pi*no_P_obs);   
                
                angle22(n) = asin( P_obs2(tt,n,2) / (no_P_obs * sin(acos( (P_obs2(tt,n,3)/no_P_obs))) ));
				p_i7(n) =  k_number^2*1*sin(acos(squeeze(P_obs2(tt,n,3))/sqrt(squeeze(P_obs2(tt,n,1))^2 + squeeze(P_obs2(tt,n,2))^2 + squeeze(P_obs2(tt,n,3))^2 ))).*exp(-1i*k_number*no_P_obs)./(4*pi*no_P_obs);   
  
				p_i(n) = p_i7(n);
				
				%Apply Babinet's principle
				integral = sum(I_end + I_s);	
				if( in(n)==0)
					integral2 = integral*exp(1i*k_number*no_P_obs);
				else
					integral2 = integral*exp(1i*k_number*no_P_obs)+p_i(n);
				end
				ps0(n)=(p_i(n)-integral2);
		
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
		
	end

end
for tt=1:n_time
	for n=1:N_obs
		for ff=1:n_freq-1
			noise_final(tt,n,ff) = 10*log10(sum(10.^(noise_final_1(tt,:,n,ff)/10)));
			SPL_dBA(tt,n,ff) = 10*log10(sum(10.^(SPL_dBA_1(tt,:,n,ff)/10)));
			
			%if((sum(10.^(-spl(tt,:,n,ff)/10)))==0)
			%	spl_total(tt,n,ff) = 0.0;
			%else
				spl_total(tt,n,ff)= -10*log10(sum(10.^(-spl(tt,:,n,ff)/10)));
			%end
			spl_total_dBA(tt,n,ff)= -10*log10(sum(10.^(-att_dBA(tt,:,n,ff)/10))+1);
			SPL_no_shield(tt,n,ff) = 10*log10(sum(10.^(noise_no_shield(tt,:,n,ff)/10)));
			SPL_no_shield_dBA(tt,n,ff) = 10*log10(sum(10.^(noise_no_shield_dBA(tt,:,n,ff)/10)));
		end
	end
	
end

%---- Calculate OSPL and OASPL and total attenuation ----
for tt=1:n_time
	for nn=1:N_obs
		OASPL_dB(tt,nn) = 10*log10(sum(10.^(noise_final(tt,nn,:)/10)));		%OSPL [dB]
		OASPL_dBA(tt,nn) = 10*log10(sum(10.^(SPL_dBA(tt,nn,:)/10)));	%OASPL [dBA]
		OSPL_no_shield(tt,nn) = 10*log10(sum(10.^(SPL_no_shield(tt,nn,:)/10)));
		OASPL_no_shield(tt,nn) = 10*log10(sum(10.^(SPL_no_shield_dBA(tt,nn,:)/10)));
	end
end

toc

pi_arr = reshape(p_i,size(X,1),size(X,2)).';

return;
%Plot contour and save video
for tt=1
	nn=1;
	for(k=1:Ny)
	for l=1:Nx
	OASPL_grid(k,l) = noise_final(tt,nn,1)-SPL_no_shield(tt,nn,1);
	nn = nn+1;
	end
	end
	colormap(jet(24));
	hold on
	pcolor(X,Y,OASPL_grid);
	xlabel('x(m)');
	ylabel('y(m)');
	set(gca,'color','none')
	shading interp;
	bar=colorbar;
	bar.Label.String = '\Delta SPL [dB]';
	plot(squeeze(pt(1,1,:,1)),squeeze(pt(1,1,:,2)),'k','lineWidth',1.5);
	hold on;
	plot(source_pos(1,1,1),source_pos(1,1,2),'m.','markerSize',12);
	axis([-2 2 -2 2]);
	xlabel('x [m]');
	ylabel('y [m]');
	hold off
	h = get(gcf,'Children');
	set(h,'FontSize',20)
	F(tt)=getframe(gcf);
	%print(gcf,'foo.png','-dpng','-r300');
	%print(gcf,'foo.eps','-depsc','-r300');
end
	

fp = fopen('Test1.txt','w');
fprintf(fp,'Obs_pos_x \t Obs_pos_y \t Obs_pos_z \t SPL_no_shield [dB] \t SPL_shield [dB] \t Delta SPL [dB]\n');
for nn=1:N_obs
fprintf(fp,'%f %f %f %f %f %f\n',P_obs(nn,1),P_obs(nn,2),P_obs(nn,3),noise_final(1,nn,1),SPL_no_shield(1,nn,1),noise_final(1,nn,1)-SPL_no_shield(1,nn,1));
end
fclose(fp);