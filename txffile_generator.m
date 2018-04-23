%%Run this file.

clear all; %#ok<CLALL>
close all;
clc;

%Set the address of file
[totalData1, HR_data1, ~] = myDoReadData...
    ("D:\4th sem\Physiological Signal Processing\SWELL Dataset\Data\pp1_c1.S00");

%% Downsampling
r = 1;
fs = totalData1.fs/r;
%HR_data = decimate(HR_data1, r);
HR_data = HR_data1(1 : (floor(length(HR_data1)/fs))*fs);

%% Locating and Plotting R-R peaks
% U -> Location of R-peaks in samples
% For accurate identification of peaks please select appropriate value for
% ...MinPeakHeight & MinPeakDistance
MinPeakHeight = 500;
MinPeakDistance = 1/4;
[~, U] = peakfinder(HR_data, fs, MinPeakHeight, MinPeakDistance);

%% Declarations
% u -> Location of R-peaks in seconds
% K -> Number of R-peaks
% T -> Length of data [0, T]
% N -> R-peaks vector in samples
% H -> History vector
[u, K, T, N, H] = init(HR_data, fs, U);

% syms t; %variable for computing integ_f & integ_df

%% file write
fileID = fopen('D:\4th sem\Physiological Signal Processing\SWELL Dataset\Data\rr_peaks_pp2-N.txt','w');
fprintf(fileID,'%4.4f\r\n', u');
fclose(fileID);