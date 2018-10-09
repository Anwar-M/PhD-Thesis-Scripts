clear all; clc; close all;
addpath('O:\MATLAB Signal Processing Files');
save_data = 0*1;
save_path = 'O:\PhD Thesis\DISSERTATION\chapter-4\figures\';
fontsize = 24;

size_le_marker = 20;

mic_info=h5read('O:\Global Optimization Source Localization\Array Benchmark Dallas\b7\ab7aCsmEss.h5', ...
'/MetaData/ArrayAttributes/microphonePositionsM');

XA = mic_info(1,:);
YA = mic_info(2,:);
ZA = mic_info(3,:);
plot(XA,YA,'k.','MarkerSize',size_le_marker);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [], [-1 1], ...
               [-1 1], -1:.5:1, -1:.5:1, fontsize, 'off', 'off', 1, 0, [], save_data, [save_path 'array2']);

run 'O:\Global Optimization Source Localization\Array Benchmark Dallas\b0\CsmReader.m';

XA = mic_info(1,:);
YA = mic_info(2,:);
ZA = mic_info(3,:);
plot(XA,YA,'k.','MarkerSize',size_le_marker);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [], [-1.1 1.1], ...
               [-1.1 1.1], -1:.5:1, -1:.5:1, fontsize, 'off', 'off', 1, 0, [], save_data, [save_path 'array1']);

% clear all;
load('O:\Global Optimization Source Localization\Inversion Exp\CSM_SINE_5000HZ.mat');

XA = mic_config(1,:);
YA = mic_config(2,:);
ZA = mic_config(3,:);
plot(XA,YA,'k.','MarkerSize',size_le_marker);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [], .5*[-1 1], ...
               .5*[-1 1], -.5:.25:.5, -.5:.25:.5, fontsize, 'off', 'off', 1, 0, [], save_data, [save_path 'array3']);
close;