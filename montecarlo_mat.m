clear; clc; close all;

rng('shuffle');
winSec = 5;             
seizure_win = [0 20];  
win_duration = seizure_win(2) - seizure_win(1); 

seizureData = {

};


nSeiz = size(seizureData,1);
mean_pv_baseline = nan(nSeiz,1);
mean_pv_seizure  = nan(nSeiz,1);
std_pv_baseline  = nan(nSeiz,1);
std_pv_seizure   = nan(nSeiz,1);
irr_pv_baseline  = nan(nSeiz,1);
irr_pv_seizure   = nan(nSeiz,1);
edf_numbers      = nan(nSeiz,1);

%% ==========================================================
%% MAIN LOOP
%% ==========================================================
disp('--- INICIANDO PROCESADO CON MONTE CARLO ---');
for s = 1:nSeiz
    
    onsetDT = datetime(seizureData{s,1},'InputFormat','dd/MM/yyyy HH:mm:ss');
    fileNum = seizureData{s,2};
    startDT = datetime(seizureData{s,3},'InputFormat','dd/MM/yyyy HH:mm:ss');
    
    edf_numbers(s) = fileNum;
    
    try
        %% LOAD MAT FILE
        matPath = fullfile('E:\1JSZ6', sprintf('1JSZ6_%d.mat', fileNum));
        loaded = load(matPath);
        
        fs = loaded.hdr.Fs;        % sampling frequency
        signal_raw = loaded.data;  % channels x samples
        
        signalData = double(signal_raw(2,:))';
        
        [b,a] = butter(4,8/(fs/2),'low');
        signal_filt = filtfilt(b,a,signalData);
        
        phase = unwrap(angle(hilbert(signal_filt)));
        instFreq = (fs/(2*pi))*[diff(phase); 0];
        instFreq(instFreq < 0) = NaN;
        instFreq(instFreq > 40) = NaN; 
        
        t = startDT + seconds((0:length(instFreq)-1)/fs);
        rel_time = seconds(t - onsetDT);
        total_len = length(t);
        
     
        winSamples = round(winSec * fs);
        s_std = movmean(instFreq, winSamples, 'omitnan');
        
        seizMask = rel_time >= seizure_win(1) & rel_time <= seizure_win(2);
        if sum(seizMask) > 0
            seiz_data = s_std(seizMask);
            mean_pv_seizure(s) = mean(seiz_data, 'omitnan');
            std_pv_seizure(s)  = std(seiz_data, 'omitnan');
            irr_pv_seizure(s)  = std(seiz_data, 'omitnan') / (abs(mean(seiz_data, 'omitnan')) + 0.01);
        end
        
        win_samples_total = round(win_duration * fs);
        num_windows_to_sample = 100;
        
        temp_mean = nan(num_windows_to_sample, 1);
        temp_std  = nan(num_windows_to_sample, 1);
        temp_irr  = nan(num_windows_to_sample, 1);
        
        valid_windows = 0;
        attempts = 0;
        max_attempts = 1000; 
        
        if total_len > win_samples_total
            while valid_windows < num_windows_to_sample && attempts < max_attempts
                attempts = attempts + 1;
                idx_start = randi(total_len - win_samples_total);
                idx_end = idx_start + win_samples_total;
                
                rand_t_start = t(idx_start);
                rand_t_end = t(idx_end);
                
                % Evitar tomar datos muy cerca de la crisis
                overlap = (rand_t_start < (onsetDT + seconds(60))) && (rand_t_end > (onsetDT - seconds(60)));
                
                if ~overlap
                    valid_windows = valid_windows + 1;
                    
                    base_data = s_std(idx_start:idx_end);
                    
                    temp_mean(valid_windows) = mean(base_data, 'omitnan');
                    temp_std(valid_windows)  = std(base_data, 'omitnan');
                    temp_irr(valid_windows)  = std(base_data, 'omitnan') / (abs(mean(base_data, 'omitnan')) + 0.01);
                end
            end
            
            if valid_windows > 0
                mean_pv_baseline(s) = mean(temp_mean(1:valid_windows), 'omitnan');
                std_pv_baseline(s)  = mean(temp_std(1:valid_windows), 'omitnan');
                irr_pv_baseline(s)  = mean(temp_irr(1:valid_windows), 'omitnan');
            end
        end
        
        fprintf('Procesado .mat %d | Fs: %.0f | Ventanas MC: %d\n', fileNum, fs, valid_windows)
        
    catch ME
        fprintf('Error en .mat %d: %s\n', fileNum, ME.message)
    end
end

valid_idx = ~isnan(mean_pv_seizure) & ~isnan(mean_pv_baseline);
N_valid = sum(valid_idx);

figure('Color', 'w', 'Position', [100 100 1100 400], 'Name', 'Resultados Baseline vs Seizure');

groups = [zeros(N_valid, 1); ones(N_valid, 1)];
x_base = ones(N_valid, 1) + (rand(N_valid, 1) - 0.5) * 0.15;
x_seiz = 2 * ones(N_valid, 1) + (rand(N_valid, 1) - 0.5) * 0.15;

% 1. MEAN
subplot(1,3,1)
data_mean = [mean_pv_baseline(valid_idx); mean_pv_seizure(valid_idx)];
boxplot(data_mean, groups, 'Labels', {'MC Baseline', 'Seizure'}, 'Colors', 'br', 'Symbol', '');
hold on;
scatter(x_base, mean_pv_baseline(valid_idx), 35, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_seiz, mean_pv_seizure(valid_idx), 35, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
title('Mean Freq (Hz)', 'FontSize', 12); 
ylabel('Value', 'FontName', 'Times', 'FontSize', 12, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 11);

% 2. STD
subplot(1,3,2)
data_std = [std_pv_baseline(valid_idx); std_pv_seizure(valid_idx)];
boxplot(data_std, groups, 'Labels', {'MC Baseline', 'Seizure'}, 'Colors', 'br', 'Symbol', '');
hold on;
scatter(x_base, std_pv_baseline(valid_idx), 35, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_seiz, std_pv_seizure(valid_idx), 35, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
title('Std Freq (Hz)', 'FontSize', 12); 
grid on; set(gca, 'FontSize', 11);

% 3. IRR 
subplot(1,3,3)
data_irr = [irr_pv_baseline(valid_idx); irr_pv_seizure(valid_idx)];
boxplot(data_irr, groups, 'Labels', {'MC Baseline', 'Seizure'}, 'Colors', 'br', 'Symbol', '');
hold on;
scatter(x_base, irr_pv_baseline(valid_idx), 35, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_seiz, irr_pv_seizure(valid_idx), 35, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
title('Irregularity (CV)', 'FontSize', 12); 
grid on; set(gca, 'FontSize', 11);
