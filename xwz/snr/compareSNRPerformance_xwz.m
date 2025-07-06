function [] = compareSNRPerformance()

%% 频率估计算法对比仿真
clear all; close all; clc;

%% 参数设置
fs = 10240;                  % 采样频率(Hz)
N = 1024;                   % 采样点数
f_useful = 505;             % 有用信号频率(Hz)
SNR_range = -10:2:30;       % 信噪比变化范围（dB)
num_runs = 50;              % 每个参数点的重复试验次数

%% ===== 第三部分：估计频率随信噪比的变化 =====
%% 初始化结果存储
est_fft = zeros(length(SNR_range), num_runs);
est_quadratic = zeros(length(SNR_range), num_runs);
est_quinn = zeros(length(SNR_range), num_runs);
est_rife = zeros(length(SNR_range), num_runs);

%% 主循环 - 改变信噪比
for i = 1:length(SNR_range)
    SNR_dB = SNR_range(i);
    
    % 对每个SNR进行多次试验以获取统计特性
    for run = 1:num_runs
        % 生成时间序列
        t = (0:N-1)/fs;
        
        % 生成有用信号
        phase_useful = 2*pi*rand;
        x_useful = A_useful * sin(2*pi*f_useful*t + phase_useful);
        
        % 生成噪声
        signal_power = mean(x_useful.^2);
        noise_power = signal_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power) * randn(size(t));
        
        % 叠加信号
        x = x_useful + noise;
        
        % 应用各种频率估计算法
        est_fft(i, run) = fft_peak_estimate(x, fs);
        est_quadratic(i, run) = quadratic_estimate(x, fs);
        est_quinn(i, run) = quinn_estimate(x, fs);
        est_rife(i, run) = rife_estimate(x, fs);
    end
end

%% 计算平均估计值
mean_est_fft = mean(est_fft, 2);
mean_est_quadratic = mean(est_quadratic, 2);
mean_est_quinn = mean(est_quinn, 2);
mean_est_rife = mean(est_rife, 2);

%% 绘制估计频率随SNR变化曲线
figure('Position', [100, 100, 800, 600]);
plot(SNR_range, mean_est_fft, 'b-o', ...
     SNR_range, mean_est_quadratic, 'r-s', ...
     SNR_range, mean_est_quinn, 'g-d', ...
     SNR_range, mean_est_rife, 'm-^', 'LineWidth', 1.5);
hold on;
plot(SNR_range, f_useful*ones(size(SNR_range)), 'k--', 'LineWidth', 2); % 真实频率
grid on;
legend('FFT直接估计', '二次多项式插值', 'Quinn A&M', 'Rife插值', '真实频率');
xlabel('信噪比 (dB)');
ylabel('估计频率 (Hz)');
title('不同频率估计算法的估计值随信噪比变化');
set(gca, 'FontSize', 12);

end