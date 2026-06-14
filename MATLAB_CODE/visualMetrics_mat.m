clear; clc; close all;

seizureData = {};


fs = 207/0.99985; 
winLen = 5;         
winStep = 2;       
fc = 8; [b,a] = butter(4, fc/(fs/2), 'low');

for s_idx = 1:size(seizureData,1)
    onsetDT = datetime(seizureData{s_idx,1}, 'InputFormat','dd/MM/yyyy HH:mm:ss');
    startTime = datetime(seizureData{s_idx,3}, 'InputFormat','dd/MM/yyyy HH:mm:ss');
    edfNum = seizureData{s_idx,2};
    
    try
        matPath = fullfile('P3', sprintf('P3_%d.mat', edfNum));
        loaded = load(matPath);
        data = loaded.data;
    
        T3 = double(data.signal(:));              % EEG signal
        t  = startTime + seconds(data.time(:));   % absolute time
    
        T3_filt = filtfilt(b,a,T3);
    
        instFreq = (fs/(2*pi)) * diff(unwrap(angle(hilbert(T3_filt))));
        tFreq = t(1:end-1);
    
        viewIdx = tFreq >= (onsetDT - seconds(60)) & ...
                  tFreq <= (onsetDT + seconds(30));
    
        t_view = tFreq(viewIdx);
        freq_view = instFreq(viewIdx);
    
        samplesWin = round(winLen * fs);
        samplesStep = round(winStep * fs);
    
        v_ratio = []; t_v = []; m_val = []; s_val = [];
    
        for i = 1:samplesStep:(length(freq_view)-samplesWin)
    
            chunk = freq_view(i:i+samplesWin);
    
            m = mean(chunk,'omitnan');
            s = std(chunk,'omitnan');
    
            m_val(end+1) = m;
            s_val(end+1) = s;
            v_ratio(end+1) = s/m;
    
            t_v(end+1) = seconds(t_view(i) - onsetDT);
        end
    
        figure('Name', sprintf('MAT %d', edfNum), ...
               'Color', 'w', 'Position', [100 100 800 600]);
    
        subplot(2,1,1);
        plot(t_v, m_val, 'LineWidth', 1.5); hold on;
        plot(t_v, s_val, 'LineWidth', 1.5);
        xline(0, '--r', 'Onset');
        ylabel('Hz');
        title(sprintf('Componentes Físicos (MAT %d)', edfNum));
        legend('Media (m)','Std Dev (s)');
        grid on;
    
        subplot(2,1,2);
        plot(t_v, v_ratio, 'k', 'LineWidth', 2);
        xline(0, '--r', 'Onset');
        ylabel('Irregularidad (s/m)');
        xlabel('Tiempo relativo al Onset (s)');
        title('Evolución de la Irregularidad Relativa');
        grid on;
    
    catch ME
        fprintf('Error en MAT %d: %s\n', edfNum, ME.message);
    end
end