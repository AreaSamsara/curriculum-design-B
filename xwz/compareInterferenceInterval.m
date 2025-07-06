function [] = compareInterferenceInterval()

%% 频率估计算法对比仿真
clear all; close all; clc;

%% 参数设置
fs = 10240;                  % 采样频率(Hz)
N = 1024;                   % 采样点数
f_useful = 505;             % 有用信号频率(Hz)
SNR_dB = 20;                % 信噪比(dB)
delta_f_range = -10:0.5:10; % 干扰与目标频率差范围(Hz)
num_runs = 50;              % 每个参数点的重复试验次数

%% ===== 第一部分：频率估计误差随干扰与目标相对位置的变化 =====
%% 初始化结果存储
errors_fft = zeros(length(delta_f_range), num_runs);
errors_quadratic = zeros(length(delta_f_range), num_runs);
errors_quinn = zeros(length(delta_f_range), num_runs);
errors_rife = zeros(length(delta_f_range), num_runs);

%% 主循环 - 改变干扰信号频率
for i = 1:length(delta_f_range)
    delta_f = delta_f_range(i);
    f_interf = f_useful + delta_f;  % 干扰信号频率
    
    % 对每个频率差进行多次试验以获取统计特性
    for run = 1:num_runs
        % 生成时间序列
        t = (0:N-1)/fs;
        
        % 生成有用信号
        A_useful = 1;
        phase_useful = 2*pi*rand;
        x_useful = A_useful * sin(2*pi*f_useful*t + phase_useful);
        
        % 生成干扰信号
        A_interf = 0.5;  % 干扰信号幅度
        phase_interf = 2*pi*rand;
        x_interf = A_interf * sin(2*pi*f_interf*t + phase_interf);
        
        % 生成噪声
        signal_power = mean(x_useful.^2);
        noise_power = signal_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power) * randn(size(t));
        
        % 叠加信号
        x = x_useful + x_interf + noise;
        
        % 应用各种频率估计算法
        f_est_fft = fft_peak_estimate(x, fs);
        f_est_quadratic = quadratic_estimate(x, fs);
        f_est_quinn = estimate_quinn(x, fs);
        f_est_rife = rife_estimate(x, fs);
        
        % 计算估计误差
        errors_fft(i, run) = f_est_fft - f_useful;
        errors_quadratic(i, run) = f_est_quadratic - f_useful;
        errors_quinn(i, run) = f_est_quinn - f_useful;
        errors_rife(i, run) = f_est_rife - f_useful;
    end
end

%% 计算平均误差
mean_errors_fft = mean(abs(errors_fft), 2);
mean_errors_quadratic = mean(abs(errors_quadratic), 2);
mean_errors_quinn = mean(abs(errors_quinn), 2);
mean_errors_rife = mean(abs(errors_rife), 2);

%% 绘制频率估计误差随相对位置变化曲线
figure('Position', [100, 100, 800, 600]);
plot(delta_f_range, mean_errors_fft, 'b-o', ...
     delta_f_range, mean_errors_quadratic, 'r-s', ...
     delta_f_range, mean_errors_quinn, 'g-d', ...
     delta_f_range, mean_errors_rife, 'm-^', 'LineWidth', 1.5);
grid on;
legend('FFT直接估计', '二次多项式插值', 'Quinn A&M', 'Rife插值');
xlabel('干扰与目标频率差 (Hz)');
ylabel('频率估计平均误差 (Hz)');
title('不同频率估计算法的误差随干扰相对位置变化');
set(gca, 'FontSize', 12);

end