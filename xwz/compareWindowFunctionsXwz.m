function [] = compareWindowFunctionsXwz()
%% Rife插值算法在不同加窗下的性能比较（500-510Hz均匀分布）
clear; clc; close all;

% 参数设置
fs = 1280;                  % 采样频率(Hz)
N = 128;                    % 样本点数
f_start = 500;              % 起始频率(Hz)
f_end = 510;                % 结束频率(Hz)
resolution = 0.5;           % 分辨率(Hz)
SNR_dB = 30;                % 信噪比
num_trials_per_freq = 100;  % 每个频率点的试验次数

% 定义测试频率点（10Hz间隔）
test_freqs = f_start:resolution:f_end;
num_freqs = length(test_freqs);

% 定义测试的窗函数
windows = {
    '矩形窗 (无窗)',   @(N) ones(1, N);
    '汉宁窗',         @hann;
    '汉明窗',         @hamming;
    '布莱克曼窗',     @blackman;
    '平顶窗',         @flattopwin;
    '凯泽窗(β=3)',    @(N) kaiser(N, 3)';
    '切比雪夫窗',     @(N) chebwin(N, 60)'
};
num_windows = size(windows, 1);

% 初始化结果存储
mse_results = zeros(num_windows, num_freqs);
bias_results = zeros(num_windows, num_freqs);

% 主测试循环（按频率点遍历）
for f_idx = 1:num_freqs
    f0 = test_freqs(f_idx);
    
    for w = 1:num_windows
        window_name = windows{w, 1};
        window_func = windows{w, 2};
        
        errors = zeros(num_trials_per_freq, 1);
        
        for trial = 1:num_trials_per_freq
            % 生成信号
            n = 0:N-1;
            signal = exp(1j * 2 * pi * f0/fs * n);
            
            % 加窗
            win = window_func(N);
            if iscolumn(win), win = win'; end  % 确保是行向量
            signal_windowed = signal .* win;
            
            % 添加噪声
            noise_power = 10^(-SNR_dB/10);
            noise = sqrt(noise_power/2) * (randn(1, N) + 1j*randn(1, N));
            x = signal_windowed + noise;
            
            % Rife插值估计
            f_est = rife_interp_windowed(x, fs, window_func);
            errors(trial) = abs(f_est - f0);
        end
        
        % 记录统计结果
        mse_results(w, f_idx) = mean(errors.^2);
        bias_results(w, f_idx) = mean(errors);
    end
    
    fprintf('已完成频率 %.1fHz 的测试\n', f0);
end

%% 结果可视化
figure('Position', [100, 100, 1400, 800]);

% 1. 各窗函数在不同频率下的MSE
subplot(2, 2, 1);
hold on;
colors = lines(num_windows);
for w = 1:num_windows
    plot(test_freqs, mse_results(w, :), 'o-', ...
        'Color', colors(w,:), 'LineWidth', 1.5, 'MarkerFaceColor', colors(w,:));
end
title('不同窗函数的MSE随频率变化');
xlabel('信号频率 (Hz)');
ylabel('均方误差 (Hz²)');
legend(windows(:,1), 'Location', 'NorthEast');
grid on;
xticks(test_freqs);

% 2. 各窗函数在不同频率下的平均误差
subplot(2, 2, 2);
hold on;
for w = 1:num_windows
    plot(test_freqs, bias_results(w, :), 's-', ...
        'Color', colors(w,:), 'LineWidth', 1.5, 'MarkerFaceColor', colors(w,:));
end
title('不同窗函数的平均误差随频率变化');
xlabel('信号频率 (Hz)');
ylabel('平均绝对误差 (Hz)');
legend(windows(:,1), 'Location', 'NorthEast');
grid on;
xticks(test_freqs);

% 3. 各窗函数性能对比（箱线图）
% subplot(2, 2, 3);
% all_errors = [];
% groups = [];
% for w = 1:num_windows
%     all_errors = [all_errors; bias_results(w, :)'];
%     groups = [groups; w*ones(num_freqs, 1)];
% end
% boxplot(all_errors, groups, 'Labels', windows(:,1));
% title('各窗函数误差分布对比');
% ylabel('绝对误差 (Hz)');
% set(gca, 'XTickLabelRotation', 45);
% grid on;

% 4. 最优窗函数选择（每个频率点）
[~, best_window_idx] = min(bias_results);
subplot(2, 2, 4);
histogram(best_window_idx, 'BinEdges', 0.5:1:num_windows+0.5);
title('各窗函数作为最优选择的频率点数量');
xlabel('窗函数索引');
ylabel('最优次数');
set(gca, 'XTick', 1:num_windows, 'XTickLabel', windows(:,1));
xtickangle(45);
grid on;

%% 性能统计表格
fprintf('\n===== 各窗函数在500-510Hz范围内的平均性能 =====\n');
fprintf('窗函数名称      | 平均MSE(Hz²) | 平均误差(Hz) | 误差标准差(Hz)\n');
fprintf('----------------|--------------|--------------|--------------\n');

for w = 1:num_windows
    avg_mse = mean(mse_results(w, :));
    avg_bias = mean(bias_results(w, :));
    std_bias = std(bias_results(w, :));
    
    fprintf('%-15s | %12.4e | %12.6f | %12.6f\n', ...
            windows{w,1}, avg_mse, avg_bias, std_bias);
end
end

%% 加窗Rife插值函数
function f_est = rife_interp_windowed(x, fs, window_func)
    N = length(x);
    
    % 加窗
    win = window_func(N);
    if iscolumn(win), win = win'; end
    x_windowed = x .* win;
    
    % 计算FFT
    X = fft(x_windowed);
    [~, k_max] = max(abs(X));
    
    % Rife插值
    k_prev = max(1, k_max-1);
    k_next = min(N, k_max+1);
    
    X_prev = abs(X(k_prev));
    X_curr = abs(X(k_max));
    X_next = abs(X(k_next));
    
    % 改进的Rife公式（防止除零）
    denominator = 4*X_curr - 2*(X_prev + X_next);
    if abs(denominator) < eps
        delta = 0;
    else
        delta = (X_next - X_prev) / denominator;
    end
    
    f_est = (k_max - 1 + delta) * fs / N;
    
    % 处理频率越界
    if f_est > fs/2
        f_est = f_est - fs;
    elseif f_est < -fs/2
        f_est = f_est + fs;
    end
end