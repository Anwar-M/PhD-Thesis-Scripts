function SPL_at_Array_Shielding_PropWing
close all;
save_data = 1;
mainresultpath = 'O:\PhD Thesis\RESULTS\CH6\PROPWING';

addpath('O:\MATLAB Signal Processing Files');

maindatapath1 = 'O:\V-Tunnel 12-12 Wing-Prop\Data MAT\';
maindatapath2 = 'O:\V-Tunnel 22-05-2018 Wing-Prop Wind\Data MAT\';
maindatapath3 = 'O:\V-Tunnel 20-11 Wing-Source\Data MAT\';
load([maindatapath1 '\..\mic_poses_optim.mat']);

fileshields = {[maindatapath1 'Meas2.mat']; [maindatapath1 'Meas4.mat']; [maindatapath1 'Meas6.mat'];...
               [maindatapath2 'Meas2.mat']; [maindatapath2 'Meas4.mat']; [maindatapath2 'Meas6.mat'];...
               [maindatapath3 'Meas01.mat']};
filesources = {[maindatapath1 'Meas1.mat']; [maindatapath1 'Meas3.mat']; [maindatapath1 'Meas5.mat'];...
               [maindatapath2 'Meas1.mat']; [maindatapath2 'Meas3.mat']; [maindatapath2 'Meas5.mat'];...
               [maindatapath3 'Meas06.mat']};
dBranges = [-4 3; -4 3; -4 3; -4 3; -4 3; -4 3; -4 3];
titlestrings = {'RPM = 4400 [min$^{-1}$]'; 'RPM = 7000 [min$^{-1}$]'; 'RPM = 7600 [min$^{-1}$]'; ...
                'RPM = 4400 [min$^{-1}$]'; 'RPM = 7000 [min$^{-1}$]';'RPM = 7600 [min$^{-1}$]'; 'Omni source'};
filenom = {'4400noflow';'7000noflow';'7600noflow';'4400flow';'7000flow';'7600flow';'omnisource'};

for I = 1:3
fileshield = fileshields{I};
filesource = filesources{I};

dBrange = dBranges(I,:);

Fs = 50e3; % sample frequency

    fprintf('\tReading shielded data...\n');
    load(fileshield);

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
    load(filesource);

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


    fl = 100;
    fu = 6300;

If = find(freq>fl);
SPLbands = 20*log10(sqrt(8/3)*sqrt(sum(psdx(If,:)*Fs/N_b/4,1))/2e-5);
% staccer(mic_poses, SPLbands, 0, [], ...
%     '$L_{\mathrm{p,exp}}$ [dB]', save_data, [mainresultpath '\' fileshield(1:end-4) '_' num2str(fc) 'Hz3rd']);
% close;

Ifns = find(freqns>fl);
SPLbandsns = 20*log10(sqrt(8/3)*sqrt(sum(psdxns(Ifns,:)*Fs/N_b/4,1))/2e-5);
% staccer(mic_poses, SPLbandsns, 0, [], ...
%     '$L_{\mathrm{p,exp}}$ [dB]', save_data, [mainresultpath '\' filesource(1:end-4) '_' num2str(fc) 'Hz3rdns']);
% close;
% 
staccer(mic_poses, SPLbands-SPLbandsns, dBrange, titlestrings{I}, ...
    '$\Delta L_{\mathrm{p}}$ [dB]', save_data, [mainresultpath '\' filenom{I} 'scatter'], I);

end


fprintf([reverseStr, '\tDone!\n']);
           
function staccer(mic_poses, dB, dBrange, strtit, ctit, save_data, file_path, I)
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
    elseif dBdiff < 7
        getickt = 1;
    else
        getickt = 2;
    end
end
if strcmp(ctit,'$\delta_{\mathrm{m}}$ [dB]') mindB = 0; end
hold on
text(mic_poses(1,Imin)-0.04, mic_poses(2,Imin)-.07, num2str(round(10*dB(Imin))/10));
text(mic_poses(1,Imax)-0.04, mic_poses(2,Imax)-.07, num2str(round(10*dB(Imax))/10));

xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
placement_x2 = 35*(xlim(2)-xlim(1))/48 + xlim(1);
placement_y2 = 112*(ylim(2)-ylim(1))/120 + ylim(1);

if I < 4
    text(placement_x2, placement_y2, ['No flow'], ...
        'Color','k','FontSize', 16,'Interpreter','LaTex', 'BackgroundColor', [1 1 1]);
elseif I < 7
    text(placement_x2, placement_y2, ['flow'], ...
        'Color','k','FontSize', 16,'Interpreter','LaTex', 'BackgroundColor', [1 1 1]);
end

hold off

c = colorbar;
h = colormap('jet');
h([2:2:24 63:-2:40],:) = [];
colormap(h);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', strtit, [-1 1], ...
               [-1 1], -1:.5:1, -1:.5:1, 16, 'on', 'on', 1, ...
               [1 mindB:getickt:round(maxdB)], ctit, save_data, file_path);