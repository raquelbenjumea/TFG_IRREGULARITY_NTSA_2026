clear; clc; close all;


seizureData = {
};

winSec = 5;

t_before = 60;
t_after  = 50;

nSeiz = size(seizureData,1);


for s = 1:nSeiz

    onsetDT = datetime(seizureData{s,1},'InputFormat','dd/MM/yyyy HH:mm:ss');
    fileNum = seizureData{s,2};
    startDT = datetime(seizureData{s,3},'InputFormat','dd/MM/yyyy HH:mm:ss');

    try
        matPath = fullfile('E:\1JSZ6', sprintf('1JSZ6_%d.mat', fileNum));
        loaded = load(matPath);

        fs = loaded.hdr.Fs;          % sampling frequency
        signal = loaded.data;        % 2 x N

        T3 = double(signal(2,:))';   % channel 2
        T3 = T3(:);                  % ensure column

        winSamples = round(winSec * fs);

        [b,a] = butter(4, 5/(fs/2), 'low');
        T3_filt = filtfilt(b,a,T3);

        instFreq = (fs/(2*pi))*[diff(unwrap(angle(hilbert(T3_filt))));0];

        t = startDT + seconds((0:length(instFreq)-1)/fs);

        s_std  = movstd(instFreq, winSamples);
        s_mean = movmean(instFreq, winSamples);

        s_dyn = s_mean;%./(abs(s_mean)+0.01);
        s_dyn(s_dyn > 10) = 10;

        idxWindow = t >= (onsetDT-seconds(t_before)) & ...
                    t <= (onsetDT+seconds(t_after));

        t_win = t(idxWindow);
        s_win = s_dyn(idxWindow);
        rel_time = seconds(t_win - onsetDT);

        if s == 1
            t_axis = -t_before:1/fs:t_after;
            irregularity_matrix = nan(nSeiz,length(t_axis));
        end

        s_interp = interp1(rel_time,s_win,t_axis,'linear',nan);

        irregularity_matrix(s,:) = s_interp;

        fprintf('Processed MAT file %d\n',fileNum)

    catch ME
        fprintf('Error in file %d: %s\n',fileNum, ME.message)
    end

end


figure('Color','w','Position',[200 200 900 500])

imagesc(t_axis,1:nSeiz,irregularity_matrix)
set(gca,'YDir','normal')

xlabel('Time relative to seizure onset (s)')
ylabel('Seizure number')
title('M aligned to seizure onset')

colorbar
colormap(jet)

hold on
xline(0,'w','LineWidth',2)


mean_irregularity = mean(irregularity_matrix,1,'omitnan');

figure
plot(t_axis,mean_irregularity,'LineWidth',3)

xlabel('Time relative to seizure onset (s)')
ylabel('M')
title('Average M across seizures')

grid on
xline(0,'r','LineWidth',2)