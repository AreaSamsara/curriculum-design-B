function [] = compareInterferenceAmtipleXwz()

%% 频率估计算法对比仿真
clear all; close all; clc;

%% 参数设置
fs = 10240;                  % 采样频率(Hz)
N = 1024;                   % 采样点数
f_useful = 505;             % 有用信号频率(Hz)
SNR_dB = 20;                % 信噪比(dB)
sir_range = 0:2:10;         % 信干比范围(dB)
num_runs = 50;              % 每个参数点的重复试验次数

%% ===== 第二部分：估计频率随信干比的变化 =====
%% 初始化结果存储
estimates_fft = zeros(length(sir_range), num_runs);
estimates_quadratic = zeros(length(sir_range), num_runs);
estimates_quinn = zeros(length(sir_range), num_runs);
estimates_rife = zeros(length(sir_range), num_runs);

%% 主循环 - 改变信干比
delta_f = 5;  % 固定频率差为5Hz
for i = 1:length(sir_range)
    sir_dB = sir_range(i);
    f_interf = f_useful + delta_f;  % 干扰信号频率
    
    % 计算幅度比
    sir_linear = 10^(sir_dB/10);
    A_useful = 1;
    A_interf = A_useful / sqrt(sir_linear);  % 调整干扰信号幅度以实现目标SIR
    
    % 对每个SIR进行多次试验以获取统计特性
    for run = 1:num_runs
        % 生成时间序列
        t = (0:N-1)/fs;
        
        % 生成有用信号
        phase_useful = 2*pi*rand;
        x_useful = A_useful * sin(2*pi*f_useful*t + phase_useful);
        
        % 生成干扰信号
        phase_interf = 2*pi*rand;
        x_interf = A_interf * sin(2*pi*f_interf*t + phase_interf);
        
        % 生成噪声
        signal_power = mean(x_useful.^2);
        noise_power = signal_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power) * randn(size(t));
        
        % 叠加信号
        x = x_useful + x_interf + noise;
        
        % 应用各种频率估计算法
        estimates_fft(i, run) = fft_peak_estimate(x, fs);
        estimates_quadratic(i, run) = quadratic_estimate(x, fs);
        estimates_quinn(i, run) = quinn_estimate(x, fs);
        estimates_rife(i, run) = rife_estimate(x, fs);
    end
end

%% 计算平均估计值
mean_estimates_fft = mean(estimates_fft, 2);
mean_estimates_quadratic = mean(estimates_quadratic, 2);
mean_estimates_quinn = mean(estimates_quinn, 2);
mean_estimates_rife = mean(estimates_rife, 2);

%% 绘制估计频率随SIR变化曲线
figure('Position', [100, 100, 800, 600]);
plot(sir_range, mean_estimates_fft, 'b-o', ...
     sir_range, mean_estimates_quadratic, 'r-s', ...
     sir_range, mean_estimates_quinn, 'g-d', ...
     sir_range, mean_estimates_rife, 'm-^', 'LineWidth', 1.5);
hold on;
plot(sir_range, f_useful*ones(size(sir_range)), 'k--', 'LineWidth', 2); % 真实频率
grid on;
legend('FFT直接估计', '二次多项式插值', 'Quinn A&M', 'Rife插值', '真实频率');
xlabel('信干比 (dB)');
ylabel('估计频率 (Hz)');
title('不同频率估计算法的估计值随信干比变化');
set(gca, 'FontSize', 12);

end