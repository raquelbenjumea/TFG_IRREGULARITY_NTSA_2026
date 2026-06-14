clc; clearvars; close all;

set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultLineLineWidth', 1.5);


rng('shuffle');
baseFolder = 'E:\';
fs = 207/0.99985;
winSecs    = 20;
winSamples = round(winSecs * fs);
baselineSecs    = 20;
baselineSamples = round(baselineSecs * fs);

bands = {
    'Delta 1-4 Hz',       [1  4];
    'Theta 4-8 Hz',       [4  8];
    'Alpha 8-12 Hz',      [8 12];
    'Beta 12-30 Hz',      [12 30];
    'Low Gamma 30-45 Hz', [30 45];
    'No Filter 1-46 Hz',  [1 46];
};
nBands = size(bands,1);
metricNames = {'Mean freq. (Hz)', 'STD (Hz)', 'V (std/mean)'};
chanNames   = {'CH1', 'CH2'};
nChans      = 2;

colCH1       = '#2563a8'; % Azul (Canal 1)
colCH2       = '#b45309'; % Naranja (Canal 2)
colBL_Seiz   = '#7c3aed'; % Morado  (Baseline vs Seizure)
colPreS_Seiz = '#b91c1c'; % Rojo    (Pre-Seizure vs Seizure)
colPreB_BL   = '#0284c7'; % Azul Claro (Pre-Baseline vs Baseline)
colPreB_PreS = '#059669'; % Verde   (Pre-Baseline vs Pre-Seizure)


seizureData_P2 = {

};

patients = {
    'P2', fullfile(baseFolder,'P2'), seizureData_P2;
};
totalPacientes = size(patients,1);

for p = 1:totalPacientes
    patientName = patients{p,1};
    dataFolder  = patients{p,2};
    seizureData = patients{p,3};
    nSeiz = size(seizureData,1);

    fprintf('\n======================================================\n');
    fprintf('>>> Procesando paciente: %s | Total de crisis: %d\n', patientName, nSeiz);
    fprintf('======================================================\n');

    patient_F = cell(nBands, nChans);
    patient_S = cell(nBands, nChans);
    patient_M = cell(nBands, nChans);

    statsReport = cell(nBands * 3 * nChans, 7);
    rowIdx = 1;

    for bIdx = 1:nBands
        bandName  = bands{bIdx,1};
        bandRange = bands{bIdx,2};
        fprintf('\n  --> Banda de Frecuencia: %s\n', bandName);

        [b_filt, a_filt] = butter(4, bandRange/(fs/2), 'bandpass');

        mat_F1 = nan(nSeiz,4); mat_S1 = nan(nSeiz,4); mat_M1 = nan(nSeiz,4);
        mat_F2 = nan(nSeiz,4); mat_S2 = nan(nSeiz,4); mat_M2 = nan(nSeiz,4);

        for s = 1:nSeiz
            edfNum  = seizureData{s,2};
            fprintf('      [Crisis %2d/%2d] (EDF: %d)... ', s, nSeiz, edfNum);
            tic;

            onsetDT = datetime(seizureData{s,1},'InputFormat','dd/MM/yyyy HH:mm:ss');
            startDT = datetime(seizureData{s,3},'InputFormat','dd/MM/yyyy HH:mm:ss');
            rel_onset = seconds(onsetDT - startDT);
            edfPath = fullfile(dataFolder, sprintf('%s_%d.edf', patientName, edfNum));

            if ~exist(edfPath,'file')
                fprintf('NO encontrado.\n'); continue;
            end

            tbl = edfread(edfPath);
            CH1 = double(vertcat(tbl{:,1}{:}));
            CH2 = double(vertcat(tbl{:,2}{:}));

            [iFreq1, t_sig] = getInstFreq(CH1, b_filt, a_filt, fs, bandRange);
            [iFreq2, ~]     = getInstFreq(CH2, b_filt, a_filt, fs, bandRange);

            % Col 1: Baseline (Interictal Late)
            mat_F1(s,1) = getMean(iFreq1, t_sig, rel_onset - 2*winSecs, rel_onset - winSecs, 0);
            mat_S1(s,1) = getStd (iFreq1, t_sig, rel_onset - 2*winSecs, rel_onset - winSecs, 0);
            mat_F2(s,1) = getMean(iFreq2, t_sig, rel_onset - 2*winSecs, rel_onset - winSecs, 0);
            mat_S2(s,1) = getStd (iFreq2, t_sig, rel_onset - 2*winSecs, rel_onset - winSecs, 0);

            % Col 2: Pre-Baseline (Interictal Early)
            mat_F1(s,2) = getMean(iFreq1, t_sig, rel_onset - 3*winSecs, rel_onset - 2*winSecs, 0);
            mat_S1(s,2) = getStd (iFreq1, t_sig, rel_onset - 3*winSecs, rel_onset - 2*winSecs, 0);
            mat_F2(s,2) = getMean(iFreq2, t_sig, rel_onset - 3*winSecs, rel_onset - 2*winSecs, 0);
            mat_S2(s,2) = getStd (iFreq2, t_sig, rel_onset - 3*winSecs, rel_onset - 2*winSecs, 0);

            % Col 3: Pre-Seizure (Ictal Early)
            mat_F1(s,3) = getMean(iFreq1, t_sig, rel_onset - winSecs, rel_onset, winSamples);
            mat_S1(s,3) = getStd (iFreq1, t_sig, rel_onset - winSecs, rel_onset, winSamples);
            mat_F2(s,3) = getMean(iFreq2, t_sig, rel_onset - winSecs, rel_onset, winSamples);
            mat_S2(s,3) = getStd (iFreq2, t_sig, rel_onset - winSecs, rel_onset, winSamples);

            % Col 4: Seizure (Ictal Late)
            mat_F1(s,4) = getMean(iFreq1, t_sig, rel_onset, rel_onset + 15, 15*fs);
            mat_S1(s,4) = getStd (iFreq1, t_sig, rel_onset, rel_onset + 15, 15*fs);
            mat_F2(s,4) = getMean(iFreq2, t_sig, rel_onset, rel_onset + 15, 15*fs);
            mat_S2(s,4) = getStd (iFreq2, t_sig, rel_onset, rel_onset + 15, 15*fs);

            mat_M1(s,:) = mat_S1(s,:) ./ abs(mat_F1(s,:));
            mat_M2(s,:) = mat_S2(s,:) ./ abs(mat_F2(s,:));

            fprintf('Ok (%.2f seg)\n', toc);
        end

        patient_F{bIdx,1} = mat_F1; patient_F{bIdx,2} = mat_F2;
        patient_S{bIdx,1} = mat_S1; patient_S{bIdx,2} = mat_S2;
        patient_M{bIdx,1} = mat_M1; patient_M{bIdx,2} = mat_M2;

    
        allMetrics = { {mat_F1, mat_F2}, {mat_S1, mat_S2}, {mat_M1, mat_M2} };

        for m = 1:3
            for chIdx = 1:nChans
                currData = allMetrics{m}{chIdx};
                validRows = ~any(isnan(currData), 2);
                cleanData = currData(validRows, :);
                nValid    = sum(validRows);

                if nValid > 2
                    N = nValid;

                  
                    Late_Inter  = cleanData(:,1);   % BL
                    Early_Inter = cleanData(:,2);   % PreB
                    Late_Ictal  = cleanData(:,4);   % Seiz
                    Early_Ictal = cleanData(:,3);   % PreS

                    Early = [Early_Inter; Early_Ictal];
                    Late  = [Late_Inter;  Late_Ictal];
                    Condition = categorical([repmat({'Interictal'}, N, 1); ...
                                            repmat({'Ictal'},      N, 1)]);
                    SubjectID = categorical(repmat((1:N)', 2, 1));

                    t_anova = table(SubjectID, Condition, Early, Late);
                    Time    = table(categorical({'Early'; 'Late'}), ...
                                   'VariableNames', {'Time'});

                    rm = fitrm(t_anova, 'Early-Late ~ Condition', ...
                               'WithinDesign', Time);

                    n_comparisons = 4;

                    % Between-subject post-hoc (Condition × Time)
                    mc_Between = multcompare(rm, 'Condition', 'By', 'Time');
                    idx_Late   = find(mc_Between.Time == 'Late',  1);
                    idx_Early  = find(mc_Between.Time == 'Early', 1);
                    p_BL_Seiz_raw   = mc_Between.pValue(idx_Late);   % BL vs Seiz
                    p_PreB_PreS_raw = mc_Between.pValue(idx_Early);  % PreB vs PreS

                    % Within-subject post-hoc (Time × Condition)
                    mc_Within  = multcompare(rm, 'Time', 'By', 'Condition');
                    idx_Inter  = find(mc_Within.Condition == 'Interictal', 1);
                    idx_Ictal  = find(mc_Within.Condition == 'Ictal',      1);
                    p_PreB_BL_raw   = mc_Within.pValue(idx_Inter);   % PreB vs BL
                    p_PreS_Seiz_raw = mc_Within.pValue(idx_Ictal);   % PreS vs Seiz

                    % Corrección de Bonferroni
                    p_BL_Seiz   = min(p_BL_Seiz_raw   * n_comparisons, 1.0);
                    p_PreB_PreS = min(p_PreB_PreS_raw  * n_comparisons, 1.0);
                    p_PreB_BL   = min(p_PreB_BL_raw    * n_comparisons, 1.0);
                    p_PreS_Seiz = min(p_PreS_Seiz_raw  * n_comparisons, 1.0);
                else
                    p_BL_Seiz   = NaN; p_PreS_Seiz = NaN;
                    p_PreB_BL   = NaN; p_PreB_PreS = NaN;
                end

                statsReport(rowIdx, :) = { ...
                    chanNames{chIdx}, bandName, metricNames{m}, ...
                    p_BL_Seiz, p_PreS_Seiz, p_PreB_BL, p_PreB_PreS };
                rowIdx = rowIdx + 1;

            end % chIdx
        end % m
    end % bIdx


    T_Stats = cell2table(statsReport, 'VariableNames', ...
        {'Channel', 'Band', 'Metric', ...
         'p_Betw_BL_Seiz', 'p_With_PreS_Seiz', ...
         'p_With_PreB_BL', 'p_Betw_PreB_PreS'});

    fprintf('\n==================================================================================================\n');
    fprintf('TABLA DE RESULTADOS ESTADÍSTICOS (BONFERRONI) - PACIENTE: %s\n', patientName);
    fprintf('==================================================================================================\n');
    disp(T_Stats);

    excelPath = fullfile(baseFolder, ...
        sprintf('Estadisticas_ANOVA_Bonferroni_%s.xlsx', patientName));
    writetable(T_Stats, excelPath);
    fprintf('\n--> Estadísticas exportadas a: %s\n', excelPath);

    condLabels = {'BL','Pre-B','Pre-S','Seiz'};
    allResults = {patient_F, patient_S, patient_M};

    for chPlot = 1:nChans
        fig = figure('Color','w','Position',[100 100 1100 650], ...
                     'Name', sprintf('Grid - %s - %s', patientName, chanNames{chPlot}));
        t_plot = tiledlayout(3, nBands, 'TileSpacing','compact','Padding','compact');

        
        chFilter = strcmp(T_Stats.Channel, chanNames{chPlot});
        T_ch = T_Stats(chFilter, :);  % subtabla solo para este canal

        for mIdx = 1:3
            for bIdx = 1:nBands
                nexttile; hold on;

                dat = allResults{mIdx}{bIdx, chPlot};
                validRows = ~any(isnan(dat), 2);
                d_clean   = dat(validRows, :);

                colActive = colCH1;
                if chPlot == 2, colActive = colCH2; end

                if ~isempty(d_clean)
                    m_vals = mean(d_clean, 1);

                    plot(1:4, m_vals, '-o', 'Color', colActive, ...
                         'LineWidth', 2.5, 'MarkerFaceColor', colActive, 'MarkerSize', 5);

                    xticks(1:4); xlim([0.5 4.5]);
                    set(gca, 'FontSize', 10, 'LineWidth', 1); grid on;

                    yMin = min(m_vals); yMax = max(m_vals);
                    yRange = yMax - yMin;
                    if yRange == 0, yRange = 0.1; end
                    ylim([yMin - 0.2*yRange, yMax + 0.85*yRange]);

                    % Recuperar p-values de T_ch para esta banda y métrica
                    rowMask = strcmp(T_ch.Band, bands{bIdx,1}) & ...
                              strcmp(T_ch.Metric, metricNames{mIdx});
                    if any(rowMask)
                        p_BL_Seiz   = T_ch.p_Betw_BL_Seiz(rowMask);
                        p_PreS_Seiz = T_ch.p_With_PreS_Seiz(rowMask);
                        p_PreB_BL   = T_ch.p_With_PreB_BL(rowMask);
                        p_PreB_PreS = T_ch.p_Betw_PreB_PreS(rowMask);
                    else
                        p_BL_Seiz = NaN; p_PreS_Seiz = NaN;
                        p_PreB_BL = NaN; p_PreB_PreS = NaN;
                    end

                    currentYOffset = 0.15;

                    % 1. Within: Pre-Baseline vs Baseline (cols 2 vs 1)
                    if ~isnan(p_PreB_BL) && p_PreB_BL < 0.05
                        sigStr = getSigStr(p_PreB_BL);
                        yPos = yMax + currentYOffset*yRange;
                        drawBracket(gca, 1, 2, yPos, yRange, colPreB_BL, sigStr);
                        currentYOffset = currentYOffset + 0.17;
                    end

                    % 2. Between: Pre-Baseline vs Pre-Seizure (cols 2 vs 3)
                    if ~isnan(p_PreB_PreS) && p_PreB_PreS < 0.05
                        sigStr = getSigStr(p_PreB_PreS);
                        yPos = yMax + currentYOffset*yRange;
                        drawBracket(gca, 2, 3, yPos, yRange, colPreB_PreS, sigStr);
                        currentYOffset = currentYOffset + 0.17;
                    end

                    % 3. Within: Pre-Seizure vs Seizure (cols 3 vs 4)
                    if ~isnan(p_PreS_Seiz) && p_PreS_Seiz < 0.05
                        sigStr = getSigStr(p_PreS_Seiz);
                        yPos = yMax + currentYOffset*yRange;
                        drawBracket(gca, 3, 4, yPos, yRange, colPreS_Seiz, sigStr);
                        currentYOffset = currentYOffset + 0.17;
                    end

                    % 4. Between: Baseline vs Seizure (cols 1 vs 4)
                    if ~isnan(p_BL_Seiz) && p_BL_Seiz < 0.05
                        sigStr = getSigStr(p_BL_Seiz);
                        yPos = yMax + currentYOffset*yRange;
                        drawBracket(gca, 1, 4, yPos, yRange, colBL_Seiz, sigStr);
                    end
                end

                if mIdx == 3
                    xticklabels(condLabels);
                else
                    xticklabels({});
                end
                if mIdx == 1
                    title(bands{bIdx,1}, 'FontSize', 12, 'FontWeight', 'bold');
                end
                if bIdx == 1
                    ylabel(metricNames{mIdx}, 'FontWeight','bold','FontSize',11);
                end
                hold off;
            end % bIdx
        end % mIdx

        title(t_plot, sprintf('%s — %s', patientName, chanNames{chPlot}), ...
              'FontSize', 13, 'FontWeight', 'bold');

        %% Exportación por canal
        pngPath = fullfile(baseFolder, ...
            sprintf('GRID_Bonferroni_%s_%s.png', patientName, chanNames{chPlot}));
        exportgraphics(fig, pngPath, 'Resolution', 600);
        pdfPath = fullfile(baseFolder, ...
            sprintf('GRID_Bonferroni_%s_%s.pdf', patientName, chanNames{chPlot}));
        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        fprintf('--> Guardado %s y PDF en: %s\n', chanNames{chPlot}, baseFolder);

    end % chPlot
end % p

disp('==========================================================');
disp('¡PROCESO COMPLETADO CON ÉXITO PARA TODOS LOS PACIENTES!');
disp('==========================================================');


%% ==========================================================
%% FUNCIONES AUXILIARES
%% ==========================================================
function s = getSigStr(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'ns';
    end
end

function drawBracket(ax, x1, x2, yPos, yRange, col, sigStr)
    tickH = 0.06 * yRange;
    xMid  = (x1 + x2) / 2;
    plot(ax, [x1 x2], [yPos yPos],      'Color', col, 'LineWidth', 1.5);
    plot(ax, [x1 x1], [yPos yPos-tickH],'Color', col, 'LineWidth', 1.5);
    plot(ax, [x2 x2], [yPos yPos-tickH],'Color', col, 'LineWidth', 1.5);
    text(ax, xMid, yPos + 0.02*yRange, sigStr, ...
         'HorizontalAlignment','center','VerticalAlignment','bottom', ...
         'Color', col, 'FontSize', 12, 'FontWeight','bold');
end

function [instFreq, t] = getInstFreq(signal, b, a, fs, bandRange)
    sig_filt = filter(b, a, signal);
    phase    = unwrap(angle(hilbert(sig_filt)));
    instFreq = (fs/(2*pi)) * [0; diff(phase)];
    instFreq(instFreq < 0 | instFreq > bandRange(2)) = NaN;
    t = (0:length(instFreq)-1) / fs;
end

function m = getMean(iFreq, t, tStart, tEnd, minSamples)
    mask = t >= tStart & t <= tEnd;
    seg  = iFreq(mask); seg = seg(~isnan(seg));
    if length(seg) > 0.5*minSamples || (minSamples == 0 && ~isempty(seg))
        m = mean(seg);
    else
        m = NaN;
    end
end

function s = getStd(iFreq, t, tStart, tEnd, minSamples)
    mask = t >= tStart & t <= tEnd;
    seg  = iFreq(mask); seg = seg(~isnan(seg));
    if length(seg) > 0.5*minSamples || (minSamples == 0 && ~isempty(seg))
        s = std(seg);
    else
        s = NaN;
    end
end

