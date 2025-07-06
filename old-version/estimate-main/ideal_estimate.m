% 无噪声无干扰条件下，估计单频正弦信号频率
% 输入：
% f_seq：信号的频率
% delta_f：频率分辨率
% N：采样点数
% window_type：
function [f_estimate,bias_error] = ideal_estimate(f_seq ,delta_f ,N ,window_type ,estimate_type)
bias_error = zeros(1,length(f_seq));    %预分配
f_estimate = zeros(1,length(f_seq));    %预分配
for n = 1:length(f_seq)
    t = (0:N-1)/(N*delta_f);    % 生成时间序列
    A = 1;                      % 信号幅度
    x = A*cos(2*pi*f_seq(n)*t);          % 信号采样序列
   %% 加窗
   switch(window_type)
       case('none')
            windows = ones(1,length(x));   %不加窗，即乘上全1矩阵
       case('hanning')
            windows = hann(length(x))';   %汉宁窗
       case('hamming')
            windows = hamming(length(x))';    %汉明窗
       case('blackman')        
            windows = blackman(length(x))';   %布莱克曼窗
       case('blackmanharris')
            windows = blackmanharris(length(x))';  %布莱克曼-哈里斯窗
   end
   signal_window = x.*windows;        %信号乘上窗函数
   %% 估计
    switch(estimate_type)
        case('direct')
            f_estimate(n) = direct_estimate(signal_window ,delta_f ,N); %直接估计
        case('quinn&AM')
            f_estimate(n) = quinn_AM(signal_window,N,N*delta_f);    %quinn-A&M算法估计
        case('jacobsen')
            f_estimate(n) = Jacobsen_Interpolation(N*delta_f,signal_window,window_type);    %jacobsen插值估计
        case('segment FFT')
            f_estimate(n) = Segment_FFT(N*delta_f,signal_window,pi/10);      %分段FFT算法估计
    end
    bias_error(n) = f_estimate(n) - f_seq(n);      % 计算bias
end
