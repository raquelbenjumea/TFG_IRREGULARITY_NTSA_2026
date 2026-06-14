clear; clc; close all
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');

rootPath    = 'E:\'; % Ruta base donde están las carpetas de pacientes
numPatients = 3;


pIDs{1}       = 'P1';
edfNums(1)    = ;
onsetStrs{1}  = ;
startStrs{1}  = ;
fs_values(1)  = 207 / 0.99985; 

pIDs{2}       = 'P2';
edfNums(2)    = ;          
onsetStrs{2}  = ; 
startStrs{2}  = ; 
fs_values(2)  = 207 / 0.99985;        

pIDs{3}       = 'P3';
edfNums(3)    = ;          
onsetStrs{3}  = ; 
startStrs{3}  = ; 
fs_values(3)  = 207 / 0.99985;        

t_before    = 600; % Segundos antes del onset
t_after     = 600; % Segundos después del onset
freqLimits  = [0.5 34]; % Límites de frecuencia para CWT

disp('Procesando y generando figura combinada (A, B, C) sin barras de color...');

fig = figure('Units', 'centimeters', 'Position', [2 2 24 22], 'Color', 'w');
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for p = 1:numPatients
    currentPID   = pIDs{p};
    currentEDF   = edfNums(p);
    currentOnset = onsetStrs{p};
    currentStart = startStrs{p};
    currentFS    = fs_values(p);
    patientPath  = fullfile(rootPath, currentPID);
    
    fprintf('\n--- Procesando %s (%d/%d) ---\n', currentPID, p, numPatients);
    
    try
        onsetDT = datetime(currentOnset, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
        startDT = datetime(currentStart, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
    catch
        error('Error en formato de fecha/hora para %s.', currentPID);
    end
    
    edfPath = fullfile(patientPath, sprintf('%s_%d.edf', currentPID, currentEDF));
    if ~exist(edfPath, 'file')
        warning('No se encuentra el archivo EDF para %s. Saltando paciente.', currentPID);
        continue; 
    end
    
    fprintf('   Leyendo archivo EDF...\n');
    tbl = edfread(edfPath);
    X = cellfun(@double, tbl{:,:}, 'UniformOutput', false);
    X = cell2mat(X);
    x = X(:, 1); 
    
    clear tbl X; 
    
    t_full   = startDT + seconds((0:length(x)-1)/currentFS);
    rel_time = seconds(t_full - onsetDT); 
    
    idx_win = rel_time >= -t_before & rel_time <= t_after;
    x_win   = x(idx_win);
    t_win   = rel_time(idx_win);
    
    clear x t_full rel_time;
    
    if length(x_win) < currentFS * 10
        warning('El archivo de %s termina antes de capturar la crisis completa. Saltando.', currentPID);
        continue;
    end
    
    x_win(isnan(x_win)) = 0; 
    
    fprintf('   Calculando Espectrograma CWT...\n');
    [cfs, f_cwt] = cwt(double(x_win), 'amor', currentFS, ...
                       'FrequencyLimits', freqLimits, ...
                       'VoicesPerOctave', 12);
    s_power = abs(cfs).^2;
    
    step = 5; 
    t_plot = t_win(1:step:end);
    s_plot = s_power(:, 1:step:end);
    
    fprintf('   Renderizando panel...\n');
    ax = nexttile;
    
    [TT, FF] = meshgrid(t_plot, f_cwt);
    surf(ax, TT, FF, zeros(size(s_plot)), s_plot, 'EdgeColor', 'none');
    shading(ax, 'flat'); 
    axis(ax, 'xy'); 
    view(ax, 2);
    
    set(ax, 'YScale', 'log', ...
            'YLim', freqLimits, ...
            'YTick', [1 4 8 12 16 30], ... 
            'YTickLabel', {'1','4','8','12','16','30'}, ...
            'TickDir', 'out', 'FontSize', 12);
    
    colormap(ax, jet); 
    caxis(ax, [0 32]); 
    xline(ax, 0, 'w--', 'LineWidth', 2); 
    
    ylabel(ax, 'Freq (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    xlim(ax, [-t_before t_after]);
    
    set(ax, 'Box', 'off'); 
    grid(ax, 'on');
    
   
    panelLetter = char(64 + p); % 65 es 'A', 66 es 'B', 67 es 'C'
    text(ax, -0.06, 1.05, panelLetter, 'Units', 'normalized', ...
         'FontSize', 16, 'FontWeight', 'bold');
         
    if p == numPatients
        xlabel(ax, 'Relative time to onset (s)', 'FontSize', 12, 'FontWeight', 'bold');
    else
        xticklabels(ax, {}); 
    end
    
    clear cfs s_power TT FF t_plot s_plot x_win t_win;
    
end 

outputFilename = 'Espectrogramas_Combinados.png';
fprintf('\nGuardando imagen combinada en alta resolución: %s...\n', outputFilename);
drawnow;
exportgraphics(fig, outputFilename, 'Resolution', 300);

fprintf('\n¡Proceso completado!\n');