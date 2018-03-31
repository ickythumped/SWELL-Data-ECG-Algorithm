function [peaks, locs] = peakfinder(HR_data, fs, min_height, min_dist)
%Identifying peaks and location of spikes
%   min_dist should be multiple of 1/2
t1 = 0:1/fs:(length(HR_data)-1)/fs;
figure
plot(t1, HR_data)
hold on
findpeaks(HR_data, 'MinPeakHeight', min_height, 'MinPeakDistance', min_dist*fs)
[peaks, locs] = findpeaks(HR_data, 'MinPeakHeight', min_height,...
    'MinPeakDistance', min_dist*fs);
hold off


end

