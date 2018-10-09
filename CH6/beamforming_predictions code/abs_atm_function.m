function alpha = abs_atm_function(f, T, H)

% This program used the program: factor_eta.m
%
% abs_atm_function: same as abd_atm, but now in function form and without
% plotting.
% This program is called for in SPL.m
%
% calculates absorption coefficient in dB/m
% as a function of
% - frequency f [Hz]
% - temperature T [degrees Celsius]
% - relative humidity H [%]


if f <= 4000
    f0 = f;
else
x = [4000 4000;
     5000 4500;
     6300 5600;
     8000 7100;
    10000 9000];
f0 = interp1(x(:,1),x(:,2),f);
end

delta = sqrt(1010./f0) * ...
( 10^(log10(H) - 1.328924 + 3.179768e-2*T) )*( 10^(-2.173716e-4*T^2 + 1.7496e-6*T^3) );

eta = factor_eta(delta);
a = 10^(2.05*log10(f0/1000) + 1.1394e-3*T - 1.916984);
b = 10^(log10(f0) + 8.42994e-3*T - 2.755624);

% alpha in dB/100m
alpha1 = a + eta*b;
% in dB/m
alpha = alpha1/100;

% end of program