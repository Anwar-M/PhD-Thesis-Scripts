% Based on EPS, run immediately after! Can only be ran once due to changing
% strings.
close all;

list = dir('*.eps');
for I = 1:numel(list)
    fid = fopen(list(I).name,'rt') ;
    X = fread(fid) ;
    fclose(fid) ;
    X = char(X.') ;
    % replace string S1 with string S2
    Y = strrep(X, '%%BoundingBox:     0     0   421   316', '%%BoundingBox:     0     0   387   316') ;
    fid2 = fopen(list(I).name,'wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
end