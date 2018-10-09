function AdjustBF
load('PROP.mat');
I = 1;

fileender = {'PROP'; 'PYLPROP'; 'PYLPROPBLOW'};
titles = {'Isolated propeller';
          'Pylon propeller';
          'Pylon propeller blowing'};
      
save_data = 0;

figure;
contourf(x, z, SPLxz, (round(maxvalxz)-dynamic_range):0.25:round(maxvalxz));

if strcmp(fileender{I}, 'PROP')
    layoutPylonprop(gca,0);
else
    layoutPylonprop(gca,1);
end

[dummy, vert_ind] = max(SPLxz);
[~, hor_ind] = max(dummy);
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
Xsh = xlim(1) + (hor_ind-1)*reso;
Zsh = ylim(1) + (vert_ind(hor_ind)-1)*reso;
placement_x = 1*(xlim(2)-xlim(1))/40 + xlim(1);
placement_y = 1*(ylim(2)-ylim(1))/30 + ylim(1);

% line([xlim(1); xlim(2)], [Zsh; Zsh], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);
% line([Xsh; Xsh], [ylim(1); ylim(2)], 'Color', 'k', 'LineStyle', '--','LineWidth', 1);

% text(placement_x, placement_y, ['Max: (' num2str(Xsh) ', ' num2str(Zsh) ')'], ...
%     'Color','k','FontSize', 14,'Interpreter','LaTex', 'BackgroundColor', 0.93*[1 1 1]);

plot_settings_font(gca, '$x$ [m]', '$z$ [m]', [titles{I}], [xmin xmax], ...
    [zmin zmax], linspace(xmin, xmax, 7), linspace(zmin, zmax, 6), 16, ...
    'on', 'on', 1, [1 round(maxvalxz)-dynamic_range/2 round(maxvalxz)], ...
    '$L_{\mathrm{p}}$ [dB]', save_data, ...
    ['..\BFxz' fileender{I} '_y=' num2str(scan_plane_Y)]);

keyboard;

function layoutPylonprop(hAxes, pylon)

% Pylon-prop config
parkerpen = 8;

% Prop
Xp = [-.37 -1.0 -.4];
Xpu = [-.37 -1.0 -0.146];
Xpd = [-.37 -1.0 -0.654];

% Pylon
Xpyl_corn1 = [-.5224 -1.0 -.4];
Xpyl_corn2 = [-.5224 -1.0 -1.3];
Xpyl_corn3 = [-1.011 -1.0 -1.3];
Xpyl_corn4 = [-1.011 -1.0 -.4];

% DNW support above
Sup1 = [Xp + [.066 0 0.05]; Xp + [.066 0 -0.05]; ...
        Xp + [1.500 0 -0.05]; Xp + [1.500 0 0.05]; ...
        Xp + [1.133+0.3 0 0.05]; Xp + [1.133+0.3 0 0.05+0.333]; ...
        Xp + [1.133+2 0 0.05+0.333]; Xp + [1.133+2 0 0.05+0.333+0.133]; ...
        Xp + [1.133 0 0.05+0.466]; Xp + [1.133 0 0.05]];
    
% DNW support below
Sup2 = [Xpyl_corn3; Xpyl_corn3 + [1.167 0 0]; Xpyl_corn3 + [1.167+1.7 0 .25-.083];
        Xpyl_corn3 + [1.167+1.7 0 -.25-.083]; Xpyl_corn3 + [1.167 0 -.166]; ...
        Xpyl_corn3 + [-.333 0 -.166]; Xpyl_corn3 + [-.333-.075 0 -.12]; ...
        Xpyl_corn3 + [-.433 0 -.083]; Xpyl_corn3 + [-.333-.075 0 -.04];...
        Xpyl_corn3 + [-.333 0 0];];

hold(hAxes, 'on')
plot(hAxes, Xp(1),Xp(3),'ko','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpu(1),Xpu(3),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);
plot(hAxes, Xpd(1),Xpd(3),'kx','MarkerSize', parkerpen, 'LineWidth', 1.5);
patch([Xpu(1); Xpd(1)], [Xpu(3); Xpd(3)], 'k', 'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
if pylon
patch([Xpyl_corn1(1); Xpyl_corn2(1); Xpyl_corn3(1); Xpyl_corn4(1); Xpyl_corn1(1)], ...
      [Xpyl_corn1(3); Xpyl_corn2(3); Xpyl_corn3(3); Xpyl_corn4(3); Xpyl_corn1(3)], 'k', ...
      'Parent', hAxes, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
end
patch(Sup1(:,1), Sup1(:,3), 'k--', 'Parent', hAxes, 'FaceColor', 'none', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
patch(Sup2(:,1), Sup2(:,3), 'k--', 'Parent', hAxes, 'FaceColor', 'none', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
hold(hAxes, 'off')
