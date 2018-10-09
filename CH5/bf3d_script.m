function bf3d_script(A, x, y, z, xs, ys, zs, dynamic_range, mic_config, doMax, doSource, xsource)
% doMax = 1;

% dynamic_range = 12*1;

SPL = 20*log10( sqrt(real(A)) / 2e-5 );
% SPL = 20*log10( sqrt(real(A)) / (4*pi*1*2e-5) );
maxdB = ceil(max(SPL(:)));
mindB = maxdB - dynamic_range;

xsize = length(x);
ysize = length(y);
zsize = length(z);
% SPL = permute(SPL,[2 3 1]);

X = repmat(ones(ysize,1)*x, 1, 1, zsize);
Y = repmat(y'*ones(1,xsize), 1, 1, zsize);
Z = permute(repmat(ones(ysize,1)*z, 1, 1, xsize),[3 1 2]);


xmin = min(X(:)); 
ymin = min(Y(:)); 
zmin = min(Z(:));

xmax = max(X(:)); 
ymax = max(Y(:)); 
zmax = max(Z(:));

figure('Position', [550 200 840 630]);

colormap(jet);
h = slice(X,Y,Z,SPL, xs, ys, zs);
set(h,'EdgeColor','none',...
'FaceColor','interp',...
'FaceAlpha','interp')
alpha('color')

hold on
h2 = slice(X,Y,Z,SPL,xmax,ymin,zmin);
set(h2,'EdgeColor','none',...
'FaceColor','interp',...
'FaceAlpha','interp')
alpha('color')
hold off

%%
hold on

% Xar = [[0 0; 0 0; 0 0; 0 0]; ...
%        zeros(8,2); zeros(8,2); ...
%        [0 0.4]; [0 -0.4]; [0 0.4]; [0 -0.4]]; 
% Yar = [[-1 1; 1 1; 1 -1; -1 -1]; ...
%        linspace(-1,1,8)'*[1 1]; ones(8,1)*[-1 1]; ...
%        [1 1]; [1 1]; [-1 -1]; [-1 -1]]; 
% Zar = [[-1 -1; -1 1; 1 1; 1 -1]; ...
%        ones(8,1)*[-1 1]; linspace(-1,1,8)'*[1 1]; ...
%        [-1 -1.572]; [-1 -1.572]; [-1 -1.572]; [-1 -1.572]];
% line(Xar', Yar', Zar', 'Color', 'k', 'LineWidth', 1.5);

scatter3(mic_config(:,1),mic_config(:,2),mic_config(:,3),'k.', 'LineWidth',1.5)

% [Xc,Yc,Zc] = cylinder([0.3 0.3]);
% Xc = Xc + 1.46;
% Zc = Zc - 1.572;
% Zc(2,:) = Zc(1,:) + 0.332;
% mesh(Xc,Yc,Zc,'facealpha',0.5, 'LineWidth', 1.5);

hold off

axis([xmin xmax ymin ymax zmin zmax]);
% axis([-1 xmax -1.6 ymax zmin zmax]);

h = gca;  % Handle to currently active axes
set(h, 'YDir', 'reverse');

caxis([mindB maxdB]);

daspect([1,1,1])
axis tight
view(-38.5,16)

ylabel('$y$ [m]');
xlabel('$x$ [m]');
zlabel('$z$ [m]');

% halpha = round(dynamic_range/1.5);
% halpha = round(dynamic_range/1);
% lalpha = 64-halpha;
% alphamap([zeros(1,lalpha) 0.625*ones(1,halpha)]);

hc = colorbar('horiz');
% set(hc, 'Interpreter', 'Latex');

if doSource
    hold on
    scatter3(xsource(:,1),xsource(:,2),xsource(:,3),'k*', 'LineWidth',1.5)
%     patch(wing_distance*[1 1 1 1], [-1.23 1.23 1.23 -1.23]/2, ...
%           [-0.0970 -0.0970 0.148 0.148], [.75 .75 .75], 'FaceAlpha', 0.5);
    hold off
end

if doMax
    hold on
    [dum, xmaxArr] = max(SPL);
    [dum2, ymaxArr] = max(dum);
    [~, zmax_ind] = max(dum2);
    xmax_ind = xmaxArr(1,ymaxArr(ymax_ind),zmax_ind);
    ymax_ind = xmaxArr(ymax_ind);

    clear dum2 dum xmaxArr zmaxArr;
    
%     Xlm = [z(zmax_ind) z(zmax_ind); z(zmax_ind) z(zmax_ind); xmin xmax];
%     Ylm = [x(xmax_ind) x(xmax_ind); ymin ymax; x(xmax_ind) x(xmax_ind)];
%     Zlm = [zmin zmax; y(ymax_ind) y(ymax_ind); y(ymax_ind) y(ymax_ind)];
    Xlm = [z(xmax_ind) z(xmax_ind); z(xmax_ind) z(xmax_ind); xmin xmax];
    Ylm = [x(zmax_ind) x(zmax_ind); ymin ymax; x(zmax_ind) x(zmax_ind)];
    Zlm = [zmin zmax; y(ymax_ind) y(ymax_ind); y(ymax_ind) y(ymax_ind)];
    
    line(Xlm', Ylm', Zlm', 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
    
    hold off
end

end