clc; clearvars; close all;

matPath   = ""; 
startTime = datetime();
onsetTime = datetime();

winSecs = 10;
stepSecs = 5;

if ~exist(matPath, 'file')
    error('No se encuentra el archivo: %s', matPath);
end

loaded = load(matPath);
fs = loaded.hdr.Fs;          
signal = loaded.data;        

F3 = double(signal(1,:))';   
T3 = double(signal(2,:))';   

totalSamples = numel(F3);
t = startTime + seconds((0:totalSamples-1)/fs);

winSamples = round(winSecs * fs);
stepSamps  = round(stepSecs * fs);

f = figure('Name','EEG Signal Viewer','Units','normalized',...
    'Position',[0.1 0.1 0.8 0.8],'KeyPressFcn',@keyControl);

uimenu(f,'Label','Ir a Tiempo','Callback',@jumpToTime);

tl = tiledlayout(2,1,'TileSpacing','compact','Padding','tight');

% Plot Canal F3
ax1 = nexttile(1);
hF3 = plot(ax1, t(1:winSamples), F3(1:winSamples), 'b', 'LineWidth', 0.8);
ylim(ax1, [-250 250]);
ylabel(ax1, 'µV');
title(ax1, 'Canal F3 - EEG');
grid on;

% Plot Canal T3
ax2 = nexttile(2);
hT3 = plot(ax2, t(1:winSamples), T3(1:winSamples), 'r', 'LineWidth', 0.8);
ylim(ax2, [-250 250]);
ylabel(ax2, 'µV');
title(ax2, 'Canal T3 - EEG');
grid on;

linkaxes([ax1 ax2], 'x');


for ax = [ax1 ax2]
    xline(ax, onsetTime, '--r', 'Seizure Onset', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
end

data = struct( ...
    'fs', fs, ...
    'winSamples', winSamples, ...
    'stepSamps', stepSamps, ...
    'currentStart', 1, ...
    'totalSamples', totalSamples, ...
    't', t, 'F3', F3, 'T3', T3, ...
    'ax1', ax1, 'ax2', ax2, ...
    'hF3', hF3, 'hT3', hT3);

guidata(f, data);


function keyControl(src, event)
    data = guidata(src);
    newStart = data.currentStart;
    
    switch event.Key
        case 'rightarrow'
            newStart = data.currentStart + data.stepSamps;
        case 'leftarrow'
            newStart = data.currentStart - data.stepSamps;
        case 'j'
            jumpToTime(src, []);
            return;
        otherwise
            return;
    end
    
    newStart = max(1, min(newStart, data.totalSamples - data.winSamples + 1));
    data.currentStart = newStart;
    guidata(src, data);
    refreshPlots(src);
end

function jumpToTime(src, ~)
    f = ancestor(src, 'figure');
    data = guidata(f);
    prompt = {'Ingrese tiempo (dd/MM/yyyy HH:mm:ss):'};
    dlgTitle = 'Saltar a Timestamp';
    defaultAns = {datestr(data.t(data.currentStart), 'dd/mm/yyyy HH:MM:SS')};
    answer = inputdlg(prompt, dlgTitle, [1 50], defaultAns);
    
    if isempty(answer), return; end
    
    try
        targetTime = datetime(answer{1}, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
        [~, idx] = min(abs(data.t - targetTime));
        data.currentStart = max(1, min(idx, data.totalSamples - data.winSamples + 1));
        guidata(f, data);
        refreshPlots(f);
    catch
        errordlg('Formato de fecha incorrecto. Use: dd/MM/yyyy HH:mm:ss');
    end
end

function refreshPlots(f)
    data = guidata(f);
    idxRange = data.currentStart : data.currentStart + data.winSamples - 1;
    
    tWin = data.t(idxRange);
    set(data.hF3, 'XData', tWin, 'YData', data.F3(idxRange));
    set(data.hT3, 'XData', tWin, 'YData', data.T3(idxRange));
    
    % Actualizar límites del eje X
    xlim(data.ax1, [tWin(1) tWin(end)]);
end