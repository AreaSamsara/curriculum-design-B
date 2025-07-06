% 交互界面demo【如果退不出去了按ctrl+C】
% 信号频率——situation.f_seq
% 干扰频率——situation.Interfere
% 噪声类型——situation.noise_type
% 迭代次数——situation.Iteration
% 信噪比/信干比——situation.SNR_dB
% 选择窗函数，直接修改windows_type
% 选择估计方法，直接修改estimate_method

addpath(genpath(pwd));
clear;close all;
old_menu();