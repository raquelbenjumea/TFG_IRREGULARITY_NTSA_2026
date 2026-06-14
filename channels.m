clc; clearvars; close all;

set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');

baseFolder  = 'E:\';
fs          = 207 / 0.99985;

colorCH1 = [0.00 0.28 0.55]; 
colorCH2 = [0.75 0.20 0.10]; 


patientName_1 = 'P2';
dataFolder_1  = fullfile(baseFolder, patientName_1);

winSec      = 20;
winSamples  = round(winSec * fs);

t_before_1 = 400;
t_after_1  = 400;
t_axis_1   = -t_before_1:1/fs:t_after_1;
Nt_1       = length(t_axis_1);

t_base_start = -200;
t_base_end   = -50;
idx_base     = (t_axis_1 >= t_base_start & t_axis_1 <= t_base_end);

metric_names = {'V (a.u)', 'S (Hz)', 'M(Hz)'};

seizureData_1 = {
};
nSeiz_1 = size(seizureData_1,1);

data_all = nan(nSeiz_1, 3, 2, Nt_1);
for s = 1:nSeiz_1
    fprintf('Procesando %s - Crisis %d/%d\n', patientName_1, s, nSeiz_1);
    edfNum  = seizureData_1{s,2};
    onsetDT = datetime(seizureData_1{s,1}, 'InputFormat','dd/MM/yyyy HH:mm:ss');
    startDT = datetime(seizureData_1{s,3}, 'InputFormat','dd/MM/yyyy HH:mm:ss');
    edfPath = fullfile(dataFolder_1, sprintf('%s_%d.edf', patientName_1, edfNum));
    
    if ~exist(edfPath, 'file')
        warning('Archivo no encontrado: %s', edfPath);
        continue;
    end
    
    tbl = edfread(edfPath);
    CH1 = double(vertcat(tbl{:,1}{:}));
    CH2 = double(vertcat(tbl{:,2}{:}));
    X = [CH1, CH2];
    nSamples = size(X,1);
    t_full = (0:nSamples-1)'/fs;
    
    onset_sec = seconds(onsetDT - startDT); 
    t_aligned = t_full - onset_sec; 
    
    for ch = 1:2
        sig = detrend(X(:,ch));
        
        phase = unwrap(angle(hilbert(sig)));
        instF = (fs/(2*pi)) * [0; diff(phase)];
        instF(instF < 0 | instF > 80) = NaN; % Limpieza
        
        MPV = movmean(instF, winSamples, 'omitnan');
        STD = movstd(instF, winSamples, 'omitnan');
        V   = STD ./ MPV;
        
        data_all(s,1,ch,:) = interp1(t_aligned, V,   t_axis_1, 'linear', NaN);
        data_all(s,2,ch,:) = interp1(t_aligned, STD, t_axis_1, 'linear', NaN);
        data_all(s,3,ch,:) = interp1(t_aligned, MPV, t_axis_1, 'linear', NaN);
    end
end

if nSeiz_1 > 0
    fig1 = figure('Color','w','Units','centimeters','Position',[2 2 16 18]);
    tl1 = tiledlayout(3,1,'TileSpacing','normal','Padding','compact');
    h_lines_1 = gobjects(2, 1);
    letras_panel_1 = {'A', 'B', 'C'}; 
    
    for m = 1:3
        ax1 = nexttile(tl1);
        hold(ax1,'on');
        
        for ch = 1:2
            if ch == 1, col = colorCH1; else, col = colorCH2; end
            
            mu = squeeze(mean(data_all(:,m,ch,:), 1, 'omitnan'));
            smooth_samples = round(10 * fs);
            mu_smooth = smoothdata(mu, 'gaussian', smooth_samples, 'omitnan');
            
            h_tmp = plot(ax1, t_axis_1, mu_smooth, 'Color', col, 'LineWidth', 1.7);
            if m == 1
                h_lines_1(ch) = h_tmp; 
            end
        end
        
        h_onset_1 = xline(ax1, 0, '--', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.1);
        h_onset_1.Alpha = 0.8;
        
        text(-0.03, 1.03, letras_panel_1{m}, 'Units', 'normalized', ...
            'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        
        ax1.Box = 'off';
        ax1.TickDir = 'out';
        ax1.LineWidth = 0.8;
        ax1.FontSize = 12;
        ax1.XLim = [-t_before_1 t_after_1];
        ax1.XTick = -400:100:400; 
        
        grid(ax1, 'on'); 
        ax1.GridAlpha = 0.08;
        ax1.XMinorGrid = 'on';
        ylabel(ax1, metric_names{m});
        
        if m == 3
            xlabel(ax1, 'Time relative to onset (s)');
        else
            ax1.XTickLabel = [];
        end
        ax1.Color = 'w';
    end
    
    lg1 = legend(nexttile(1), [h_lines_1(1), h_lines_1(2), h_onset_1], ...
        {'Distal-to-Central', 'Central-to-Proximal', 'Seizure onset'}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11, ...
        'FontWeight', 'bold', 'TextColor', [0.1 0.1 0.1], ...
        'Location', 'southoutside');
    lg1.Units = 'normalized';
    lg1.Position = [0.5 - lg1.Position(3)/2, 0.02, lg1.Position(3), lg1.Position(4)];
    lg1.ItemTokenSize = [18 18];
    
    set(fig1, 'PaperPositionMode', 'auto');
    exportName_1 = sprintf('FIG_THESIS_%s_Broadband_Estetica.pdf', patientName_1);
    exportgraphics(fig1, fullfile(baseFolder, exportName_1), 'ContentType', 'vector');
    disp(['Parte 1 terminada. Exportado: ', exportName_1]);
else
    disp('Parte 1 omitida: No hay crisis en seizureData_1.');
end


smoothSecs    = 2;
smoothSamples = round(smoothSecs * fs);
t_start_2     = -60;
t_end_2       = 60;
t_axis_2      = t_start_2:1/fs:t_end_2;
Nt_2          = length(t_axis_2);

bands = {
    'Delta',     [1 4];
    'Theta',     [4 8];
    'Alpha',     [8 12];
    'Beta',      [12 30];
    'Low gamma', [30 45];
    'Broadband', []
};

seizureData_P3 = {
};
patients = {'P3', fullfile(baseFolder,'P3'), seizureData_P3};

letras_panel_2 = {'A', 'B', 'C', 'D', 'E', 'F'};

for p = 1:size(patients,1)
    patientName_2 = patients{p,1};
    dataFolder_2  = patients{p,2};
    seizureData_2 = patients{p,3};
    nSeiz_2       = size(seizureData_2,1);
    
    fprintf('\n>>> Patient: %s | Seizures: %d\n', patientName_2, nSeiz_2);
    if nSeiz_2 == 0, continue; end
    
    fig2 = figure('Color','w','Units','centimeters');
    fig2.Position = [2 2 16 10];   
    tl2 = tiledlayout(3,2,'TileSpacing','normal','Padding','compact');
    ax_handles = gobjects(size(bands,1),1);
    h_lines_2  = gobjects(2,size(bands,1));
    
    for bIdx = 1:size(bands,1)
        bandName  = bands{bIdx,1};
        bandRange = bands{bIdx,2};
        fprintf(' -> Band: %s\n', bandName);
        
        all_M1 = nan(nSeiz_2, Nt_2);
        all_M2 = nan(nSeiz_2, Nt_2);
        
        if ~isempty(bandRange)
            [b_filt,a_filt] = butter(4,bandRange/(fs/2),'bandpass');
        end
        
        for s = 1:nSeiz_2
            edfNum  = seizureData_2{s,2};
            onsetDT = datetime(seizureData_2{s,1},'InputFormat','dd/MM/yyyy HH:mm:ss');
            startDT = datetime(seizureData_2{s,3},'InputFormat','dd/MM/yyyy HH:mm:ss');
            rel_onset = seconds(onsetDT - startDT);
            
            edfPath = fullfile(dataFolder_2,sprintf('%s_%d.edf',patientName_2,edfNum));
            if ~exist(edfPath,'file'), continue; end
            
            tbl = edfread(edfPath);
            CH1 = double(vertcat(tbl{:,1}{:}));
            CH2 = double(vertcat(tbl{:,2}{:}));
            
            if ~isempty(bandRange)
                sig1 = filtfilt(b_filt,a_filt,CH1);
                sig2 = filtfilt(b_filt,a_filt,CH2);
                cutoffFreq = bandRange(2);
            else
                sig1 = CH1;
                sig2 = CH2;
                cutoffFreq = fs/2;
            end
            
            [instF1,t_full2] = getContInstFreq(sig1,fs,cutoffFreq);
            [instF2,~]       = getContInstFreq(sig2,fs,cutoffFreq);
            
            M1 = movmean(instF1,smoothSamples,'omitnan');
            M2 = movmean(instF2,smoothSamples,'omitnan');
            
            rel_time = t_full2 - rel_onset;
            
            all_M1(s,:) = interp1(rel_time,M1,t_axis_2,'linear',NaN);
            all_M2(s,:) = interp1(rel_time,M2,t_axis_2,'linear',NaN);
        end
        
        mu_M1 = mean(all_M1,1,'omitnan');
        mu_M2 = mean(all_M2,1,'omitnan');
        
        ax2 = nexttile(tl2);
        ax_handles(bIdx) = ax2;
        hold(ax2,'on');
        
        h_lines_2(1,bIdx) = plot(ax2,t_axis_2,mu_M1,'Color',colorCH1,'LineWidth',1.7);
        h_lines_2(2,bIdx) = plot(ax2,t_axis_2,mu_M2,'Color',colorCH2,'LineWidth',1.7);
        
        h_onset_2 = xline(ax2,0,'--','Color',[0.85 0.1 0.1],'LineWidth',1.1);
        h_onset_2.Alpha = 0.8;
        
        text(ax2, -0.07, 1.07, letras_panel_2{bIdx}, 'Units', 'normalized', ...
            'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        
        ax2.Box = 'off';
        ax2.TickDir = 'out';
        ax2.LineWidth = 0.8;
        ax2.FontSize = 12;
        ax2.XLim = [-60 60];
        ax2.XTick = -60:20:60;
        
        grid(ax2,'on');
        ax2.GridAlpha = 0.08;
        ax2.XMinorGrid = 'on';
        
        ylabel(ax2,sprintf('%s (Hz)',bandName));
        
        if bIdx < 5
            ax2.XTickLabel = [];
        else
            xlabel(ax2,'Time relative to onset (s)');
        end
        
        ax2.Color = 'w';
        hold(ax2,'off');
    end
    
    % Leyenda Parte 2
    lg2 = legend(ax_handles(end), ...
        [h_lines_2(1,end), h_lines_2(2,end), h_onset_2], ...
        {'Distal-to-Central','Central-to-Proximal','Seizure onset'}, ...
        'Orientation','horizontal', 'Box','off', 'FontSize',11, ...             
        'FontWeight','bold', 'TextColor',[0.1 0.1 0.1], ... 
        'Location','southoutside');
    lg2.Units = 'normalized';
    lg2.Position = [0.5-lg2.Position(3)/2, 0.02, lg2.Position(3), lg2.Position(4)];
    lg2.ItemTokenSize = [18 18];  
    
    set(fig2,'PaperPositionMode','auto');
    exportName_2 = sprintf('FIG_THESIS_%s_multiband_M.pdf',patientName_2);
    exportgraphics(fig2, fullfile(baseFolder, exportName_2), 'ContentType','vector');
    disp(['Exportado: ', exportName_2]);
end

disp('=== SCRIPT FINALIZADO ===');

function [instFreq,t] = getContInstFreq(signal,fs,maxFreq)
    phase = unwrap(angle(hilbert(signal)));
    instFreq = (fs/(2*pi))*[0; diff(phase)];
    instFreq(instFreq < 0 | instFreq > maxFreq+10) = NaN;
    t = (0:length(instFreq)-1)/fs;
end