function SPL_at_Array_Shielding_Plate
close all;
save_data = 1;
mainresultpath = 'O:\PhD Thesis\RESULTS\CH6\SOURCEPLATE_SCATTER';

addpath('O:\MATLAB Signal Processing Files');

maindatapath = 'O:\V-Tunnel 30-10 Plate-Source\Data MAT';
load([maindatapath '\..\mic_poses_optim.mat']);

fileshields = {'Mes1.mat';'Mes1.mat';'Mes1.mat';...
               'Mes6.mat';'Mes6.mat';'Mes6.mat';'Mes1.mat';'Mes6.mat'};
filesources = {'Mes3.mat';'Mes3.mat';'Mes3.mat';...
               'Mes4.mat';'Mes4.mat';'Mes4.mat';'Mes3.mat';'Mes4.mat'};
filecomps = {'Meas1_2000.dat';'Meas1_4000.dat';'Meas1_5000.dat';...
             'Meas6_2000.dat';'Meas6_4000.dat';'Meas6_5000.dat';'Meas1_OSPL.dat';'Meas6_OSPL.dat'};
fcs = [2000; 4000; 5000; 2000; 4000; 5000];
dBranges = [0 3; 0 10; 0 10; 0 3; 0 10; 0 10; 0 4; 0 4];

for I = 1:6
fileshield = fileshields{I};
filesource = filesources{I};
filecomp = filecomps{I};
if I < 7 fc = fcs(I); end
dBrange = dBranges(I,:);

Fs = 50e3; % sample frequency

if (I==1)||(I==4)||(I==7)||(I==8)

    fprintf('\tReading shielded data...\n');
    load([maindatapath '\' fileshield]);

    [N, n_mic] = size(data);
    t = N/50e3; % sec
    t_b = 5*0.5; % sec
    N_b = t_b*Fs;

    n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
    fprintf('\tUsing %d blocks\n', n_blocks);
    df = Fs/N_b;

    psdx = zeros(floor(N_b/2), n_mic);
    freq = (0:floor(N_b/2)-1)*df;
    reverseStr = '';
    for B = 1:n_blocks
        msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        xdft = zeros(N_b, n_mic);
        for M = 1:n_mic
            % fft with hanning window to all microphones, index according to 50% overlap
            xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
                hann(N_b));
        end
        psdx = psdx + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

    end
    psdx = 4*psdx/n_blocks; % Hanning AMPLITUDE CORRECTION factor
    psdx(2:end-1,:) = 2*psdx(2:end-1,:); % Single-side

    fprintf([reverseStr, '\tAll blocks evaluated!\n']);

    fprintf('\tReading no-shielded data...\n');
    load([maindatapath '\' filesource]);

    [N, n_mic] = size(data);
    t = N/50e3; % sec
    t_b = 5*0.5; % sec
    N_b = t_b*Fs;

    n_blocks = 2*floor(t/t_b)-1; % amount of blocks, 50 % overlap
    fprintf('\tUsing %d blocks\n', n_blocks);
    df = Fs/N_b;

    psdxns = zeros(floor(N_b/2), n_mic);
    freqns = (0:floor(N_b/2)-1)*df;
    reverseStr = '';
    for B = 1:n_blocks
        msg = sprintf('\tEvaluating block %d/%d...\n', B, n_blocks);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        xdft = zeros(N_b, n_mic);
        for M = 1:n_mic
            % fft with hanning window to all microphones, index according to 50% overlap
            xdft(:,M) = fft(data( (B-1)*N_b/2 + 1 : (B+1)*N_b/2 , M) .* ...
                hann(N_b));
        end
        psdxns = psdxns + (1/(Fs*N_b)) * abs(xdft(1:floor(N_b/2),:)).^2;

    end
    psdxns = 4*psdxns/n_blocks; % Hanning AMPLITUDE CORRECTION factor
    psdxns(2:end-1,:) = 2*psdxns(2:end-1,:); % Single-side

    fprintf([reverseStr, '\tAll blocks evaluated!\n']);

end

if I < 7
    fl = fc*2^(-1/6);
    fu = fc*2^(1/6);
    stringtitles = ['$f =$ ' num2str(fc) ' [Hz]'];
else
    stringtitles = ['OSPL'];
    fl = 500;
    fu = 6300;
end

If = find((freq>fl)&(freq<fu));
SPLbands = 20*log10(sqrt(8/3)*sqrt(sum(psdx(If,:)*Fs/N_b/4,1))/2e-5);
% staccer(mic_poses, SPLbands, 0, [], ...
%     '$L_{\mathrm{p,exp}}$ [dB]', save_data, [mainresultpath '\' fileshield(1:end-4) '_' num2str(fc) 'Hz3rd']);
% close;

Ifns = find((freqns>fl)&(freqns<fu));
SPLbandsns = 20*log10(sqrt(8/3)*sqrt(sum(psdxns(Ifns,:)*Fs/N_b/4,1))/2e-5);
% staccer(mic_poses, SPLbandsns, 0, [], ...
%     '$L_{\mathrm{p,exp}}$ [dB]', save_data, [mainresultpath '\' filesource(1:end-4) '_' num2str(fc) 'Hz3rdns']);
% close;
% 
% staccer(mic_poses, SPLbands-SPLbandsns, 0, ['$f =$ ' num2str(fc) ' [Hz]'], ...
%     '$\Delta L_{\mathrm{p,exp}}$ [dB]', save_data, [mainresultpath '\' fileshield(1:end-4) 'diff_' num2str(fc) 'Hz3rd']);

fid = fopen(['K:\co\ance\ance-ShieldingComparison\Results\plate_results\comp data files\' filecomp]);
fgets(fid);
A = textscan(fid,'%d%f%f%f%f%f%f');
fclose(fid);
% staccer(mic_poses, (A{1,7}.'-(SPLbands-SPLbandsns)), 0, ...
%     ['$f =$ ' num2str(fc) ' [Hz]'], '$\Delta$ [dB]', ...
%     save_data, [mainresultpath '\' fileshield(1:end-4) 'diffdiff_' num2str(fc) 'Hz3rd_noabs']);
staccer(mic_poses, abs(A{1,7}.'-(SPLbands-SPLbandsns)), dBrange, ...
    stringtitles, '$\delta_{\mathrm{m}}$ [dB]', ...
    save_data, [mainresultpath '\' fileshield(1:end-4) 'diffdiff_' num2str(fc) 'Hz3rd']);

end

for I = 7:8
    Fs = 50e3;
    fl = 500;
    fu = 6300;
    
    fileshield = fileshields{I};
    filesource = filesources{I};
    filecomp = filecomps{I};
    if I < 7 fc = fcs(I); end
    dBrange = dBranges(I,:);
    
    load([maindatapath '\' fileshield]);
    
    [B,A] = butter(4,[fl fu]/(Fs/2)); % designs a highpass filter.

    for K = 1:size(data,2)
        data(:,K) = filtfilt(B,A,data(:,K));
        SPLdBS(K) = 20*log10(rms(data(:,K))/2e-5);
    end
    
    clear data;
    load([maindatapath '\' filesource]);
    [B,A] = butter(4,[fl fu]/(Fs/2)); % designs a highpass filter.

    for K = 1:size(data,2)
        data(:,K) = filtfilt(B,A,data(:,K));
        SPLdBNS(K) = 20*log10(rms(data(:,K))/2e-5);
    end
    
    fid = fopen(['K:\co\ance\ance-ShieldingComparison\Results\plate_results\comp data files\' filecomp]);
    fgets(fid);
    A = textscan(fid,'%d%f%f%f%f%f%f');
    fclose(fid);
    
    staccer(mic_poses, abs(A{1,7}.'-(SPLdBS-SPLdBNS)), dBrange, ...
        ['OSPL'], '$\delta_{\mathrm{m}}$ [dB]', save_data, ...
        [mainresultpath '\' fileshield(1:end-4) 'diffdiff_OSPL']);
    
end


fprintf([reverseStr, '\tDone!\n']);
           
function staccer(mic_poses, dB, dBrange, strtit, ctit, save_data, file_path)
figure;
scatter(mic_poses(1,:), mic_poses(2,:), [], dB, 'filled', 'MarkerEdgeColor', 'k');
[mindB, Imin] = min(dB);
[maxdB, Imax] = max(dB);
if numel(dBrange)==2
    mindB = dBrange(1);
    maxdB = dBrange(2);
    dBdiff = dBrange(2)-dBrange(1);
    if dBdiff < 4
        getickt = 0.5;
    elseif dBdiff < 6
        getickt = 1;
    else
        getickt = 2;
    end
end
if strcmp(ctit,'$\delta_{\mathrm{m}}$ [dB]') mindB = 0; end
hold on
text(mic_poses(1,Imin), mic_poses(2,Imin)+.06, num2str(round(10*dB(Imin))/10));
text(mic_poses(1,Imax), mic_poses(2,Imax)+.06, num2str(round(10*dB(Imax))/10));
hold off

c = colorbar;
h = colormap('jet');
h([2:2:24 63:-2:40],:) = [];
colormap(h);


plot_settings_font(gca, '$x$ [m]', '$y$ [m]', strtit, [-1 1], ...
               [-1 1], -1:.5:1, -1:.5:1, 16, 'on', 'on', 1, ...
               [1 mindB:getickt:round(maxdB)], ctit, save_data, file_path);