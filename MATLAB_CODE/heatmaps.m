clear; clc; close all;


seizureData = {
};

fs       = 207/0.99985;
t_before = 600;
t_after  = 600;

t_axis = -t_before : 1/fs : t_after;
nSeiz  = size(seizureData,1);

winSec     = 20;
winSamples = round(winSec * fs);


all_m = nan(nSeiz, length(t_axis));   % MPV
all_s = nan(nSeiz, length(t_axis));   % STD
all_v = nan(nSeiz, length(t_axis));   % s/m

for sIdx = 1:nSeiz
    onsetDT = datetime(seizureData{sIdx,1}, 'InputFormat','dd/MM/yyyy HH:mm:ss');
    edfNum  = seizureData{sIdx,2};
    startDT = datetime(seizureData{sIdx,3}, 'InputFormat','dd/MM/yyyy HH:mm:ss');
    fprintf('Seizure %d/%d (EDF %d)\n', sIdx, nSeiz, edfNum);
    try
        edfPath = fullfile('E:\P3', sprintf('P3_%d.edf', edfNum));
        tbl = edfread(edfPath);
        
        x = double(vertcat(tbl{:,1}{:}));
        N = length(x);
        t_full = startDT + seconds((0:N-1)/fs);
        idxWindow = t_full >= (onsetDT - seconds(t_before)) & ...
                    t_full <= (onsetDT + seconds(t_after));
        sig = x(idxWindow);
        rel_time = seconds(t_full(idxWindow) - onsetDT);
        
        [B,A] = butter(4, [1 45]/(fs/2), 'bandpass');
        sig_filt = filtfilt(B, A, sig);

        phase = unwrap(angle(hilbert(sig_filt)));
        instFreq = (fs/(2*pi)) * [diff(phase); 0];
        
        instFreq(instFreq < 0 | instFreq > 100) = NaN;
        
        % Ventanas centradas omitiendo NaNs
        m_raw = movmean(instFreq, winSamples, 'omitnan');
        s_raw = movstd(instFreq, winSamples, 'omitnan'); % ¡Ahora evalúa instFreq!
        v_raw = s_raw ./ m_raw;
        
        % Limpieza de valores extremos de V (Igual que en el Código 1)
        threshold_V = prctile(v_raw(~isnan(v_raw)), 99);
        if isnan(threshold_V), threshold_V = 5; end
        v_raw(v_raw > threshold_V) = NaN;

        all_m(sIdx,:) = interp1(rel_time, m_raw, t_axis, 'linear', NaN);
        all_s(sIdx,:) = interp1(rel_time, s_raw, t_axis, 'linear', NaN);
        all_v(sIdx,:) = interp1(rel_time, v_raw, t_axis, 'linear', NaN);
        
        fprintf('  OK\n');
    catch ME
        fprintf('  ERROR EDF %d: %s\n', edfNum, ME.message);
    end
end

figure('Color','w','Position',[100 100 1200 900])

t = tiledlayout(3, 1, 'TileSpacing', 'normal', 'Padding', 'compact');

metrics = {all_v, all_s, all_m};
labels_derecha = {'V (a.u.)', 'S (Hz)', 'M (Hz)'};
letras_panel = {'A', 'B', 'C'}; % Letras para cada panel

for k = 1:3
    nexttile
    
    dat = metrics{k};
    clim_lo = prctile(dat(:),5);
    clim_hi = prctile(dat(:),95);
    
    imagesc(t_axis, 1:nSeiz, dat, [clim_lo clim_hi])
    set(gca,'YDir','normal')
    colormap jet
    
    % Configurar el colorbar
    cb = colorbar;
    cb.Label.String = labels_derecha{k};
    %cb.Label.FontWeight = 'bold';
    cb.Label.FontSize = 12;
    
    hold on
    xline(0,'w-','LineWidth',2)
    
    text(-0.03, 1.03, letras_panel{k}, 'Units', 'normalized', ...
        'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    
  
    ylabel('Seizure #', 'FontSize', 12)
   
    
    if k == 3
        xlabel('Time (s)', 'FontSize', 12)
    else
        set(gca,'XTickLabel',[])
    end
end