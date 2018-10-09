% Script understanding wave propagation of monopole with convection, result
% not so good. Probably need ray tracing or something
clearvars;

c = 1;
t1 = 1; % delay 1
t2 = 1.5; % delay 2
u = [0 .85];

p1 = [0 0];

th = 0:2*pi/100:2*pi;
r = ones(1,numel(th))*t1*c;
[p2(:,1), p2(:,2)] = pol2cart(th,r);

for I = 1:size(p2,1)
    unit_vectors = p2(I,:)/norm(p2(I,:));
    r2(I) = t1*(c+dot(p2(I,:)/norm(p2(I,:)),u));
end
[p3(:,1), p3(:,2)] = pol2cart(th,r2);

hold on
plot(p1(1),p1(2),'rx');
plot(p2(:,1),p2(:,2),'b--');
plot(p3(:,1),p3(:,2),'b.-');
hold off
axis equal

%%
t = 0:0.01:1;
r = 0.01*ones(1,numel(th));
[cx,cy] = pol2cart(th, 1*ones(1,numel(th)));
[x, y] = pol2cart(th,r);
x2 = x;
y2 = y;
h = axes('XLim',[-1.2 1.2], 'YLim', [-1.2 1.2]);
axis equal;
hold on
plot(p1(1),p1(2),'rx');
for I = 1:numel(t)
    x = x + cx*0.01; 
    y = y + cy*0.01;
    y2 = y2 + (cy+2*0.5)*0.01;
    gwa = plot(h, x, y, 'r.');
    mwa = plot(h, x, y2,'b.');
    pause(0.1/2);
    delete(gwa);
    delete(mwa);
end
gwa = plot(h, x, y, 'r.');
mwa = plot(h, x, y2,'b.');
hold off