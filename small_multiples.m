clear; clc; close all;
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');

patients = {};
nPatients = size(patients, 1);

fs = 207 / 0.99985;
winSec = 20; 
winSamples = round(winSec * fs);
t_before = 600; t_after = 600;
t_axis = -t_before : 1/fs : t_after;
Nt = length(t_axis);
bands = {'Delta', [1 4]; 'Theta', [4 8]; 'Alpha', [8 12]; 'Beta', [12 30]; 'Low gamma', [30 45]; 'Broadband', [1 48]};
nBands = size(bands,1);
metric_names = {'V (a.u.)', 'S (Hz)', 'M (Hz)'}; 

all_mu = cell(nPatients, 3);
for p = 1:nPatients
    patientID = patients{p,1};
    patientPath = patients{p,2};
    seizureData = patients{p,3};
    nSeizures = size(seizureData, 1);
    
    tmpFile = fullfile(patientPath, sprintf('%s_%d.edf', patientID, seizureData{1,2}));
    if isfile(tmpFile)
        tbl = edfread(tmpFile);
        nChannels = min(2, size(tbl,2));
    else
        nChannels = 2; % Fallback
    end
    
    data = nan(nSeizures, nBands, nChannels, 3, Nt);
    
    for s = 1:nSeizures
        fprintf('Procesando %s - Convulsión %d/%d\n', patientID, s, nSeizures);
        
        edfNum = seizureData{s,2};
        edfPath = fullfile(patientPath, sprintf('%s_%d.edf', patientID, edfNum));
        if ~isfile(edfPath), continue; end
        tbl = edfread(edfPath);
        X = cell2mat(cellfun(@double, tbl{:,:}, 'UniformOutput', false));
        
        onset_sec = seconds(datetime(seizureData{s,1}, 'InputFormat','dd/MM/yyyy HH:mm:ss') - ...
                            datetime(seizureData{s,3}, 'InputFormat','dd/MM/yyyy HH:mm:ss'));
        onset = round(onset_sec * fs) + 1;
        s0 = onset - round(t_before * fs);
        s1 = onset + round(t_after * fs);
        
        if s0 < 1 || s1 > size(X,1), continue; end
        
        for ch = 1:nChannels
            for b = 1:nBands
                xf = bandfilter(X(s0:s1, ch), bands{b,2}, fs);
                ph = unwrap(angle(hilbert(xf)));
                iF = (fs/(2*pi)) * [diff(ph); 0];
                
                iF(iF < 0 | iF > 100) = NaN;
                
                MPV = movmean(iF, winSamples, 'omitnan');
                STD = movstd(iF, winSamples, 'omitnan');
                V   = STD ./ MPV;
                
                threshold_V = prctile(V(~isnan(V)), 99);
                if isnan(threshold_V), threshold_V = 5; end
                V(V > threshold_V) = NaN;
                
                n = min(length(MPV), Nt);
                
                data(s,b,ch,1,1:n) = V(1:n);
                data(s,b,ch,2,1:n) = STD(1:n);
                data(s,b,ch,3,1:n) = MPV(1:n);
            end
        end
    end
    
    for m = 1:3
        all_mu{p,m} = squeeze(nanmean(data(:,:,:,m,:), 1));
    end
end

ch_colors = [0.00 0.28 0.55;
             0.75 0.20 0.10];
             
for p = 1:nPatients
    fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 36 20], 'Name', patients{p,1});
    
    tl = tiledlayout(3, nBands, 'TileSpacing','normal', 'Padding','compact');
    
    for m = 1:3
        for b = 1:nBands
            ax = nexttile; hold(ax, 'on');
            mu_mat = squeeze(all_mu{p,m}(b, :, :));
            
            h_lines = gobjects(1,2);
            for ch = 1:size(mu_mat, 1)
                sig = mu_mat(ch, :);
                h_lines(ch) = plot(ax, t_axis, sig, '-', 'LineWidth', 1.7, 'Color', ch_colors(ch,:));
            end
            
            h_onset = xline(ax, 0, '--', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.1);
            h_onset.Alpha = 0.8;
            
            xlim(ax, [-100 100]);
            
            if m == 1
                title(ax, bands{b,1}, 'FontWeight','bold', 'FontSize', 12);
            end
            
            if b == 1
                ylabel(ax, metric_names{m}, 'FontWeight','bold', 'FontSize', 12);
            end
            
            ax.YAxis.TickLabelFormat = '%.2g';
            
            % Ticks adaptados a la nueva ventana
            xticks(ax, -100:50:100);
            
            if m == 3
                xlabel(ax, 'Time relative to onset (s)', 'FontSize', 11);
            else
                xticklabels(ax, {});
            end
            
            grid(ax, 'on'); 
            ax.GridAlpha = 0.08;
            ax.XMinorGrid = 'on';
            
            set(ax, 'FontSize', 11, 'Box','off', 'TickDir','out', 'LineWidth', 0.8);
        end
    end
    
    lg = legend([h_lines(1), h_lines(2), h_onset], ...
        {'Distal-to-Central','Central-to-Proximal','Seizure onset'}, ...
        'Orientation','horizontal', 'Box','off', ...
        'FontSize', 12, 'FontWeight','bold');
    lg.Layout.Tile = 'south';
    
    exportgraphics(fig, sprintf('Fig_Dynamics_%s_Landscape.pdf', patients{p,1}), ...
        'ContentType','vector', 'BackgroundColor','white');
end

function xf = bandfilter(x, band, fs)
    if isempty(band), xf = x; return; end
    [B,A] = butter(4, band/(fs/2), 'bandpass');
    xf = filtfilt(B,A,x);
end