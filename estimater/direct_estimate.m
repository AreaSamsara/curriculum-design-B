%% ֱ�ӹ��ƺ���
function [f_estimate] = direct_estimate(signal_window,delta_f,N)
     s_fft = abs(fft(signal_window,N)); %FFT�������߷���
    [~,index] = max(s_fft);             %����������������
    f_estimate = (index-1)*delta_f;     %�������Ƶ��
end