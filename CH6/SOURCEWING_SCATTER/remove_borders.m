close all;

list = dir('*.fig');
for I = 1:numel(list)
    open(list(I).name);
    set(gcf,'Units','centimeters');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) 12.5 pos(4)]);
    print(['.\' list(I).name(1:end-4)], '-depsc', '-r300');
    close;
end