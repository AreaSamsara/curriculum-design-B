%% 直接估计函数
function [f_estimate] = direct_estimate(signal_window,delta_f,N)
     s_fft = abs(fft(signal_window,N)); %FFT后求谱线幅度
    [~,index] = max(s_fft);             %搜索幅度最大的谱线
    f_estimate = (index-1)*delta_f;     %计算估计频率
end