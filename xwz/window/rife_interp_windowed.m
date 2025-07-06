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