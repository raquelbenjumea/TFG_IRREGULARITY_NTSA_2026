clear; clc; close all;

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

for s = 1:nSeiz
    
    onsetDT = datetime(seizureData{s,1},'InputFormat','dd/MM/yyyy HH:mm:ss');
    fileNum = seizureData{s,2};
    startDT = datetime(seizureData{s,3},'InputFormat','dd/MM/yyyy HH:mm:ss');
    
    edf_numbers(s) = fileNum;
    
    try
        matPath = fullfile('E:\1JSZ6', sprintf('1JSZ6_%d.mat', fileNum));
        loaded = load(matPath);
        
        fs = loaded.hdr.Fs;        % sampling frequency
        signal = loaded.data;      % channels x samples
        
        signalData = double(signal(1,:))';
        
        [b,a] = butter(4,8/(fs/2),'low');
        signal_filt = filtfilt(b,a,signalData);
        
        instFreq = (fs/(2*pi))*[diff(unwrap(angle(hilbert(signal_filt))));0];
        
        t = startDT + seconds((0:length(instFreq)-1)/fs);
        
        rel_time = seconds(t - onsetDT);
        
        idx_base = (rel_time >= -120) & (rel_time <= -60);
        mean_pv_baseline(s) = mean(instFreq(idx_base), 'omitnan');
        std_pv_baseline(s)  = std(instFreq(idx_base), 'omitnan');
        irr_pv_baseline(s)  = std_pv_baseline(s) ./ mean_pv_baseline(s);
        
      
        idx_seiz = (rel_time >= 0) & (rel_time <= 25);
        mean_pv_seizure(s) = mean(instFreq(idx_seiz), 'omitnan');
        % FIXED BUG HERE: changed idx_base to idx_seiz for the next two lines
        std_pv_seizure(s)  = std(instFreq(idx_seiz), 'omitnan');
        irr_pv_seizure(s)  = std_pv_seizure(s) ./ mean_pv_seizure(s);
        
        fprintf('Processed MAT %d\n',fileNum)
        
    catch ME
        fprintf('Error in MAT %d: %s\n',fileNum, ME.message)
    end
end

figure('Color', [0.96 0.96 0.96], 'Position', [100 200 1200 400]);

% --- Plot 1: Mean ---
subplot(1,3,1);
plotData_mean = [mean_pv_baseline, mean_pv_seizure]';
plot([1,2], plotData_mean, 'o-', 'LineWidth',1.5,'MarkerSize',6);
xlim([0.5,2.5]);
xticks([1,2]);
xticklabels({'Baseline (-120:-60s)','Seizure (0:+15s)'});
ylabel('Mean V');
title('Mean V (T3)','FontWeight','bold');
grid on;
set(gca,'Color','w','GridAlpha',0.4);

% --- Plot 2: Standard Deviation ---
subplot(1,3,2);
plotData_std = [std_pv_baseline, std_pv_seizure]';
plot([1,2], plotData_std, 'o-', 'LineWidth',1.5,'MarkerSize',6);
xlim([0.5,2.5]);
xticks([1,2]);
xticklabels({'Baseline (-120:-60s)','Seizure (0:+15s)'});
ylabel('STD V');
title('Standard Deviation V (T3)','FontWeight','bold');
grid on;
set(gca,'Color','w','GridAlpha',0.4);

% --- Plot 3: Irregularity ---
subplot(1,3,3);
plotData_irr = [irr_pv_baseline, irr_pv_seizure]';
plot([1,2], plotData_irr, 'o-', 'LineWidth',1.5,'MarkerSize',6);
xlim([0.5,2.5]);
xticks([1,2]);
xticklabels({'Baseline (-120:-60s)','Seizure (0:+15s)'});
ylabel('Irregularity V');
title('Irregularity V (T3)','FontWeight','bold');
grid on;
set(gca,'Color','w','GridAlpha',0.4);


disp(' ');
disp('--- RESULTS TABLE 1: MEAN ---');
MeanTable = table(edf_numbers, mean_pv_baseline, mean_pv_seizure, ...
    'VariableNames', {'File','Baseline_Mean_V','Seizure_Mean_V'});
disp(MeanTable);

disp('--- RESULTS TABLE 2: STANDARD DEVIATION ---');
StdTable = table(edf_numbers, std_pv_baseline, std_pv_seizure, ...
    'VariableNames', {'File','Baseline_Std_V','Seizure_Std_V'});
disp(StdTable);

disp('--- RESULTS TABLE 3: IRREGULARITY ---');
IrrTable = table(edf_numbers, irr_pv_baseline, irr_pv_seizure, ...
    'VariableNames', {'File','Baseline_Irr_V','Seizure_Irr_V'});
disp(IrrTable);