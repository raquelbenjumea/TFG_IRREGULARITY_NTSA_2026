clear; clc; close all;


seizureData = {
};

fs = 207/0.99985;


t_before     = 120;
t_after      = 60;
baseline_win = [-120 -60];
seizure_win  = [0 20];


winSec     = 5;
winSamples = round(winSec * fs);
[b, a]     = butter(4, 3/(fs/2), 'low');
nSeiz      = size(seizureData, 1);


metrics = table();
metrics.baselineMean = nan(nSeiz, 1);
metrics.seizureMean  = nan(nSeiz, 1);
metrics.baselineMin  = nan(nSeiz, 1);
metrics.seizureMin   = nan(nSeiz, 1);
metrics.baselineMax  = nan(nSeiz, 1);
metrics.seizureMax   = nan(nSeiz, 1);

for s = 1:nSeiz
    onsetDT = datetime(seizureData{s,1}, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
    edfNum  = seizureData{s,2};
    startDT = datetime(seizureData{s,3}, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
    
    try
        edfPath = fullfile('E:\P1', sprintf('P1_%d.edf', edfNum));
        tbl = edfread(edfPath);
        T3  = double(vertcat(tbl{:,2}{:}));
        
        T3_filt  = filtfilt(b, a, T3);
        phase    = unwrap(angle(hilbert(T3_filt)));
        instFreq = (fs/(2*pi)) * [diff(phase); 0];
        
        instFreq(instFreq < 0)  = NaN;
        instFreq(instFreq > 40) = NaN;
        
        t      = startDT + seconds((0:length(instFreq)-1) / fs);
        s_mean = movstd(instFreq, winSamples, 'omitnan')./movmean(instFreq, winSamples, 'omitnan');
        
        idxWindow = t >= (onsetDT - seconds(t_before)) & ...
                    t <= (onsetDT + seconds(t_after));
        
        t_win    = t(idxWindow);
        s_win    = s_mean(idxWindow);
        rel_time = seconds(t_win - onsetDT);
        
        baseMask = rel_time >= baseline_win(1) & rel_time <= baseline_win(2);
        seizMask = rel_time >= seizure_win(1)  & rel_time <= seizure_win(2);
        
        if sum(baseMask) < 20 || sum(seizMask) < 20
            fprintf('Skipping seizure %d (insufficient data)\n', s)
            continue
        end
        
        bVals = s_win(baseMask);
        sVals = s_win(seizMask);
        
        metrics.baselineMean(s) = mean(bVals, 'omitnan');
        metrics.seizureMean(s)  = mean(sVals, 'omitnan');
        metrics.baselineMin(s)  = min(bVals);
        metrics.seizureMin(s)   = min(sVals);
        metrics.baselineMax(s)  = max(bVals);
        metrics.seizureMax(s)   = max(sVals);
        
        fprintf('Processed EDF %d\n', edfNum)
    catch
        fprintf('Error in EDF %d\n', edfNum)
    end
end

validIdx = ~isnan(metrics.baselineMean);

mNames = {'Mean', 'Min', 'Max'};
dataStruct.Mean = [metrics.baselineMean(validIdx), metrics.seizureMean(validIdx)];
dataStruct.Min  = [metrics.baselineMin(validIdx),  metrics.seizureMin(validIdx)];
dataStruct.Max  = [metrics.baselineMax(validIdx),  metrics.seizureMax(validIdx)];

statResults = struct();

for i = 1:length(mNames)
    name = mNames{i};
    bData = dataStruct.(name)(:,1);
    sData = dataStruct.(name)(:,2);
    
 
    [hB, pShapiroB] = swtest(bData, 0.05);
    [hS, pShapiroS] = swtest(sData, 0.05);
    
    if hB == 0 && hS == 0
        [~, p_val] = ttest2(bData, sData);
        testType = 'T-Test (Param)';
    else
        p_val = ranksum(bData, sData);
        testType = 'Wilcoxon (Non-Param)';
    end
    
    sd_pooled = sqrt((std(bData)^2 + std(sData)^2) / 2);
    d_val = (mean(bData) - mean(sData)) / sd_pooled;
    
    statResults(i).Metrica = name;
    statResults(i).Test_Used = testType;
    statResults(i).P_Value = p_val;
    statResults(i).Significant = p_val < 0.05;
    statResults(i).Cohen_D = d_val;
    statResults(i).SW_P_Base = pShapiroB;
    statResults(i).SW_P_Seiz = pShapiroS;
    statResults(i).Mean_Baseline = mean(bData);
    statResults(i).Mean_Seizure  = mean(sData);
end


ResumenTable = struct2table(statResults);
disp(' ')
disp('===== ANALISIS ESTADISTICO CON SHAPIRO-WILK =====')
disp(ResumenTable)

fig = figure('Name','Resultados Estadísticos','Position', [100 100 1000 200],...
             'MenuBar','none','ToolBar','none','NumberTitle','off');
uitable('Parent', fig, 'Data', table2cell(ResumenTable), ...
        'ColumnName', ResumenTable.Properties.VariableNames, ...
        'Units', 'normalized', 'Position', [0 0 1 1], 'RowName', []);