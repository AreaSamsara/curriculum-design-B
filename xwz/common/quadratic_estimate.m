function f_est = quadratic_estimate(x, fs)
    % 二次多项式插值频率估计
    % 输入: x - 输入信号 
    % 输出: f_est - 估计频率(Hz) 
    N = length(x);

    % 窗函数处理
    win = hann(N)';
    x_win = x .* win;

    % 计算FFT
    X = fft(x_win);
    
    % 找到幅度谱最大值位置
    [~, k_max] = max(abs(X));

    k_prev = max(1, k_max-1);
    k_next = min(N, k_max+1);
    
    X_prev = abs(X(k_prev));
    X_curr = abs(X(k_max));
    X_next = abs(X(k_next));
    
    delta = (X_next - X_prev) / (4*X_curr - 2*(X_prev + X_next));
    f_est = (k_max - 1 + delta) * fs / N;
end