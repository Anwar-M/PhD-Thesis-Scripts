clear all; clc; close all;
addpath('O:\MATLAB Signal Processing Files');
save_data = 1;
save_path = 'O:\PhD Thesis\DISSERTATION\chapter-6\Figures\';
fontsize = 16;

size_le_marker = 20;
load('O:\V-Tunnel 13-12 Thrust Experiment\mic_poses_optim.mat');

XA = mic_poses(1,:);
YA = mic_poses(2,:);
ZA = mic_poses(3,:);
plot(XA,YA,'k.','MarkerSize',size_le_marker);

plot_settings_font(gca, '$x$ [m]', '$y$ [m]', [], [-1 1], ...
               [-1 1], -1:.5:1, -1:.5:1, fontsize, 'off', 'off', 1, 0, [], save_data, [save_path 'arrayoptim']);
