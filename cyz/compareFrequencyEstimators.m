function [] = compareFrequencyEstimators()

%% 不同噪声分布下频率估计算法性能分析（线性坐标）
clear all; close all; clc;

% 参数设置
fs = 20000;          % 采样频率 (Hz)
N = 1024;            % 采样点数
f0 = 505;            % 信号频率 (Hz)
use_window = true;   % 是否使用窗函数
snr_range = 0:1:25; % 信噪比范围 (dB)
num_trials = 500;    % 每个信噪比下的试验次数

% 噪声类型
noise_types = {'Gaussian', 'Uniform', 'Rayleigh', 'Exponential'};
colors = {'b', 'r', 'g', 'm'};
line_styles = {'-', '--', ':', '-.'};
algorithm_names = {'FFT直接估计', 'Rife插值', '二次多项式插值', 'DFT Phase'};

% 初始化结果存储
results_rmse = struct();
results_mean = struct();
results_std = struct();
for i = 1:length(noise_types)
    results_rmse.(noise_types{i}) = zeros(4, length(snr_range));
    results_mean.(noise_types{i}) = zeros(4, length(snr_range));
    results_std.(noise_types{i}) = zeros(4, length(snr_range));
end

% 创建信号
t = (0:N-1)/fs;
signal = cos(2*pi*f0*t);

% 主测试循环 (不同噪声类型)
for noise_idx = 1:length(noise_types)
    noise_type = noise_types{noise_idx};
    fprintf('\n测试噪声类型: %s\n', noise_type);
    
    % 不同信噪比下的测试
    for snr_idx = 1:length(snr_range)
        snr = snr_range(snr_idx);
        fprintf('  处理 SNR = %d dB...\n', snr);
        
        % 初始化误差存储
        errors = zeros(4, num_trials);
        
        % 多次试验取平均
        for trial = 1:num_trials
            % 生成带噪信号
            noisy_signal = add_noise(signal, snr, noise_type);
            
            % 使用四种算法估计频率
            f_est_fft = fft_peak_estimate(noisy_signal, fs, use_window);
            f_est_rife = rife_estimate(noisy_signal, fs, use_window);
            f_est_jacobsen = jacobsen_estimator(noisy_signal, fs, use_window);
            f_est_phase = dft_phase_estimation(noisy_signal, fs, use_window);
            
            % 计算绝对误差 (Hz)
            errors(1, trial) = abs(f_est_fft - f0);
            errors(2, trial) = abs(f_est_rife - f0);
            errors(3, trial) = abs(f_est_jacobsen - f0);
            errors(4, trial) = abs(f_est_phase - f0);
        end
        
        % 计算统计指标
        rmse = sqrt(mean(errors.^2, 2));
        mean_err = mean(errors, 2);
        std_err = std(errors, 0, 2);
        
        % 存储结果
        results_rmse.(noise_type)(:, snr_idx) = rmse;
        results_mean.(noise_type)(:, snr_idx) = mean_err;
        results_std.(noise_type)(:, snr_idx) = std_err;
    end
end

% 设置绘图参数
plot_width = 900;
plot_height = 700;
marker_size = 2;
line_width = 2;

% 绘制不同噪声类型下的RMSE
figure('Position', [100, 100, plot_width, plot_height]);
for algo_idx = 1:4
    subplot(2, 2, algo_idx);
    hold on;
    grid on;
    box on;
    
    title(sprintf('%s 算法', algorithm_names{algo_idx}), 'FontSize', 14);
    xlabel('信噪比 (dB)', 'FontSize', 12);
    ylabel('RMSE (Hz)', 'FontSize', 12);
    
    % 添加FFT分辨率参考线
%     fft_res = fs/N;
%     plot(xlim, [fft_res, fft_res], 'k:', 'LineWidth', 1.5);
%     text(max(snr_range), fft_res, sprintf(' FFT分辨率=%.1fHz', fft_res), ...
%         'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'k');
    
    % 绘制不同噪声类型的曲线
    for noise_idx = 1:length(noise_types)
        noise_type = noise_types{noise_idx};
        rmse = results_rmse.(noise_type)(algo_idx, :);
        
        plot(snr_range, rmse, ...
            'LineStyle', line_styles{noise_idx}, ...
            'Color', colors{noise_idx}, ...
            'LineWidth', line_width, ...
            'Marker', 'o', ...
            'MarkerSize', marker_size, ...
            'MarkerFaceColor', colors{noise_idx}, ...
            'DisplayName', noise_type);
    end
    
    % 设置纵轴范围
    max_rmse = max(results_rmse.(noise_types{1})(algo_idx, :));
    for i = 2:length(noise_types)
        max_rmse = max(max_rmse, max(results_rmse.(noise_types{i})(algo_idx, :)));
    end
    ylim([0, max_rmse * 1.1]);
    
    legend('show', 'Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 11);
end

% 绘制不同噪声类型下的平均误差
figure('Position', [100, 100, plot_width, plot_height]);
for algo_idx = 1:4
    subplot(2, 2, algo_idx);
    hold on;
    grid on;
    box on;
    
    title(sprintf('%s 算法', algorithm_names{algo_idx}), 'FontSize', 14);
    xlabel('信噪比 (dB)', 'FontSize', 12);
    ylabel('平均绝对误差 (Hz)', 'FontSize', 12);
    
%     % 添加FFT分辨率参考线
%     plot(xlim, [fft_res, fft_res], 'k:', 'LineWidth', 1.5);
    
    % 绘制不同噪声类型的曲线
    for noise_idx = 1:length(noise_types)
        noise_type = noise_types{noise_idx};
        mean_err = results_mean.(noise_type)(algo_idx, :);
        
        plot(snr_range, mean_err, ...
            'LineStyle', line_styles{noise_idx}, ...
            'Color', colors{noise_idx}, ...
            'LineWidth', line_width, ...
            'Marker', 's', ...
            'MarkerSize', marker_size, ...
            'MarkerFaceColor', colors{noise_idx}, ...
            'DisplayName', noise_type);
    end
    
    % 设置纵轴范围
    max_mean = max(results_mean.(noise_types{1})(algo_idx, :));
    for i = 2:length(noise_types)
        max_mean = max(max_mean, max(results_mean.(noise_types{i})(algo_idx, :)));
    end
    ylim([0, max_mean * 1.1]);
    
    legend('show', 'Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 11);
end

% 显示统计表格
fprintf('\n===== 各算法在不同噪声下的最大平均绝对误差 (Hz) =====\n');
fprintf('信噪比(dB)\tFFT Peak\tRife\t\tJacobsen\tPhase Diff\n');
for snr_idx = 1:length(snr_range)
    snr = snr_range(snr_idx);
    fprintf('%d\t\t', snr);
    for algo_idx = 1:4
        max_mean = 0;
        for noise_idx = 1:length(noise_types)
            noise_type = noise_types{noise_idx};
            mean_err = results_mean.(noise_type)(algo_idx, snr_idx);
            if mean_err > max_mean
                max_mean = mean_err;
            end
        end
        fprintf('%.2f\t\t', max_mean);
    end
    fprintf('\n');
end

end

%% 噪声添加函数
function noisy_signal = add_noise(signal, snr_db, noise_type)
    % 计算信号功率
    signal_power = mean(signal.^2);
    
    % 计算所需噪声功率
    snr_linear = 10^(snr_db/10);
    noise_power = signal_power / snr_linear;
    
    % 生成噪声
    N = length(signal);
    
    switch lower(noise_type)
        case 'gaussian'
            noise = randn(1, N); % 高斯噪声
        case 'uniform'
            noise = 2*rand(1, N) - 1; % 均匀分布噪声 [-1, 1]
        case 'rayleigh'
            noise = raylrnd(1, 1, N); % 瑞利分布噪声
        case 'exponential'
            noise = exprnd(1, 1, N); % 指数分布噪声
        otherwise
            error('未知噪声类型: %s', noise_type);
    end
    
    % 调整噪声功率
    noise = noise - mean(noise); % 去除直流分量
    current_noise_power = mean(noise.^2);
    noise = sqrt(noise_power / current_noise_power) * noise;
    
    % 添加噪声
    noisy_signal = signal + noise;
end

%% FFT峰值检索算法
function f_est = fft_peak_estimate(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
        coherent_gain = sum(win)/N;
    else
        win = ones(1, N);
        x_win = x;
        coherent_gain = 1;
    end
    
    % FFT计算
    X = fft(x_win, N);
    P = abs(X)/N;
    P = P / coherent_gain;
    
    % 单边频谱
    P_single = P(1:floor(N/2)+1);
    P_single(2:end-1) = 2*P_single(2:end-1);
    f_axis = (0:floor(N/2)) * fs / N;
    
    % 找到最大谱线
    [~, k0] = max(P_single);
    f_est = f_axis(k0);
end

%% RIFE算法
function f_est = rife_estimate(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
        coherent_gain = sum(win)/N;
    else
        win = ones(1, N);
        x_win = x;
        coherent_gain = 1;
    end
    
    % FFT计算
    X = fft(x_win, N);
    P = abs(X)/N;
    P = P / coherent_gain;
    
    % 单边频谱
    P_single = P(1:floor(N/2)+1);
    P_single(2:end-1) = 2*P_single(2:end-1);
    f_axis = (0:floor(N/2)) * fs / N;
    
    % 找到最大谱线
    [~, k0] = max(P_single);
    
    % 确定次大谱线方向
    if k0 == 1
        r = 1;
        mag_second = P_single(k0+1);
    elseif k0 == length(P_single)
        r = -1;
        mag_second = P_single(k0-1);
    else
        left_mag = P_single(k0-1);
        right_mag = P_single(k0+1);
        if right_mag > left_mag
            r = 1;
            mag_second = right_mag;
        else
            r = -1;
            mag_second = left_mag;
        end
    end
    
    % 计算幅度比
    mag_max = P_single(k0);
    alpha = mag_second / mag_max;
    
    % Rife算法核心公式 (汉宁窗专用)
    delta = (2*alpha - 1) / (alpha + 1);
    
    % 频率估计
    f_est = f_axis(k0) + r * delta * (fs / N);
end

%% 二次多项式算法
function f_est = jacobsen_estimator(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
        coherent_gain = sum(win)/N;
    else
        win = ones(1, N);
        x_win = x;
        coherent_gain = 1;
    end
    
    % FFT计算
    X = fft(x_win, N);
    P = abs(X)/N;
    P = P / coherent_gain;
    
    % 单边频谱
    P_single = P(1:floor(N/2)+1);
    P_single(2:end-1) = 2*P_single(2:end-1);
    f_axis = (0:floor(N/2)) * fs / N;
    
    % 找到最大谱线
    [~, k0] = max(P_single);
    
    % 提取峰值点及相邻两点幅值
    alpha = P_single(k0);
    
    if k0 > 1
        beta = P_single(k0-1);
    else
        beta = 0;
    end
    
    if k0 < length(P_single)
        gamma = P_single(k0+1);
    else
        gamma = 0;
    end
    
    % 二次多项式插值公式
    if k0 > 1 && k0 < length(P_single)
        delta = (gamma - beta) / (2*(2*alpha - beta - gamma));
    else
        delta = 0; % 边界情况使用基本FFT
    end
    
    % 频率估计
    f_est = f_axis(k0) + delta * (fs / N);
end

%% DFT相位差分法
function f_est = dft_phase_estimation(signal, fs, use_window)
% 输入:
%   signal - 输入实正弦波信号
%   fs     - 采样频率 (Hz)
% 输出:
%   f0_hat     - 估计频率 (Hz)

N = length(signal);
if mod(N, 2) ~= 0
    error('信号长度N必须是偶数');
end
% 1. 将信号分成两个等长子序列
half = N/2;
s1 = signal(1:half);     % 前一半序列
s2 = signal(half+1:end); % 后一半序列
% 窗函数处理
    if use_window
        win = hann(N/2)';
        x1 = signal(1:N/2) .* win;
        x2 = signal(N/2+1:end) .* win;
    else
        x1 = signal(1:N/2);
        x2 = signal(N/2+1:end);
    end

% 2. 分别计算两个子序列的DFT
S1 = fft(x1);  % N/2点DFT
S2 = fft(x2);  % N/2点DFT

% 3. 只考虑正频率部分（避免负频率干扰）
S1 = S1(1:floor(half/2)+1);
S2 = S2(1:floor(half/2)+1);

% 4. 找到DFT幅度最大谱线的位置（在正频率范围内）
[~, idx1] = max(abs(S1));  % 第一段最大幅度位置
[~, idx2] = max(abs(S2));  % 第二段最大幅度位置

% 5. 检查两段的最大峰值位置是否一致（防止噪声干扰）
if idx1 ~= idx2
    % 如果不一致，选择两段中幅度和更大的位置
    if (abs(S1(idx1)) + abs(S2(idx1))) > (abs(S1(idx2)) + abs(S2(idx2)))
        k0 = idx1;
    else
        k0 = idx2;
    end
else
    k0 = idx1;
end
k0 = k0 - 1;  % MATLAB索引从1开始，转换为从0开始的索引

% 6. 提取最大谱线处的相位
phi1 = angle(S1(k0+1));  % 第一段DFT相位
phi2 = angle(S2(k0+1));  % 第二段DFT相位

% 7. 计算相位差并解模糊（保持在[-π, π]范围内）
delta_phi = phi2 - phi1;
delta_phi = mod(delta_phi + pi, 2*pi) - pi;  % 相位解缠绕

% 8. 计算子序列频率分辨率
delta_f_sub = fs / half;  % 子序列DFT频率分辨率

% 9. 估计相对频偏δ
delta = delta_phi / (2*pi);

% 10. 估计频率
f_est = (k0 + delta) * delta_f_sub;

end