function [SPL_obs,freq_doppler,k_number,t_obs] = atm_effects(tt,ff,P_obs,source_pos,U,delta_t,c,freq,PSD,T,hum)

%==============================================================================================================
%Function doppler shifted the frequencies, apply spherical spreading and absoption

%Input  - tt , int					: time step
%		- ff , int					: frequency step
%		- P_obs, 2D array			: type (N observers,3) with the coordinates of the observers [m]
%		- source_pos, 2D array 		: type (tt,3) with the coordinates of the source in each time step [m]
%		- U , 1D array				: velocity of the aircraft and source (Ux,Uy,Uz) [m/s]
%		- delta_t , double 			: value of each time step [s]
%		- c , double 				: speed of sound [m/s]
%		- freq, array				: array with the frequencies of the source (Hz)
%		- PSD, array				: array with the PSD [Pa2/Hz]
%		- T , double				: temperature [deg celsius]
%		- hum, double 				: humidity [%]

%Output - SPL_obs , double			: final valor of SPL [dB]
%		- freq_doppler, 1D array	: doppler shifted frequencies
%		- k_number , 1D array		: wave number of each freq_doppler
%		- t_obs, double				: time at the observer
%==============================================================================================================



r = norm(P_obs-source_pos);
	
u = P_obs-source_pos;

v = U/norm(U);
	
mach = U/c;
	
t_obs = delta_t*(tt-1) + r/c; %time at the observer
			
theta = acos(dot(u,v)/r); % angle between flight path and observer
	
freq_doppler = freq;%/(1-norm(mach)*cos(theta)); %Dopler shifted frequencies

k_number = (2*pi*freq(ff))/c;
	
SPL = PSD; %Convert PSD (dB/Hz) to dB
	
absorp = abs_atm_function(freq_doppler(ff),T,hum); %Calculate absoption coefficient for each frequency
			
if(isnan(absorp)==1)
	absorp = abs_atm_function(10000,T,hum);
end

			
SPL_obs = SPL - 20*log10(r);% - absorp*r; %Transmission loss in the atmosphere
	