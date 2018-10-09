clear all; clc; close all;
addpath('O:\MATLAB Signal Processing Files');
save_data = 1;
save_path = 'O:\PhD Thesis\DISSERTATION\chapter-5\figures\';
fontsize = 24;

size_le_marker = 12;
load('O:\APIAN\APIAN Pylon Processed\Vortex Shedding\2014-09-16_11-34-38\configuration\config.txt');

XA = config(:,3);
YA = config(:,4);
% ZA = ARRAY_ROT(:,3);
plot(-XA+2,YA-1.5,'k.','MarkerSize',size_le_marker);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [], [-1.7 1.7], ...
               [-1.7 1.7], -1.5:0.75:1.5, -1.5:0.75:1.5, fontsize, 'off', 'off', 1, 0, [], save_data, [save_path 'arrayTUDLOG']);
% close;