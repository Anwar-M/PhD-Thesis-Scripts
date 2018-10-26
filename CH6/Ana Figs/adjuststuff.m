addpath('O:\MATLAB Signal Processing Files');
list = dir('*.fig');
save_data = 1;

dbranges = [-3 3; 0 12; -8 -2; -3 3; -4 2];
ctit = {'$\delta$ [dB]'; '$\Delta L_{\mathrm{p, exp}}$ [dB]'; ...
     '$\Delta L_{\mathrm{p, exp}}$ [dB]'; '$\delta$ [dB]'; ...
     '$\Delta L_{\mathrm{p, exp}}$ [dB]'};
file = {'array_thrust_corrected';'array_shield_OSPL_static_85';'array_propeller_airflow_85';'array_thrust_corrected_2';'array_shield_source'};
for I = 1:length(list)
    fh = open(list(I).name);
    set(fh,'WindowStyle','normal');
    set(fh, 'Position', [800         400         560         420]);
    h = findobj(gca);
    h(end).MarkerEdgeColor = [0, 0, 0];
    h(end).SizeData = 36;
%     set(gcf, 'PaperPosition', [87.6378  263.4449  420.0000  315.0000]);
    
    set(gcf,'Units','centimeters');
    pos = get(gcf,'Position');
    set(gcf,'Position',[400 400 12.5 pos(4)]);
    plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [], [-1 1], ...
               [-1 1], -1:.5:1, -1:.5:1, 16, 'on', 'on', 1, ...
               [1 dbranges(I,1) dbranges(I,2)], ctit{I}, 0, ['.\' file{I}]);
           
    set(gcf,'Units','centimeters');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) 12.5 pos(4)]);
    print(['.\' file{I}], '-depsc', '-r300');
   close;
end