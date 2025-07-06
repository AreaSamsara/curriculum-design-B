% �������޸��������£����Ƶ�Ƶ�����ź�Ƶ��
% ���룺
% f_seq���źŵ�Ƶ��
% delta_f��Ƶ�ʷֱ���
% N����������
% window_type��
function [f_estimate,bias_error] = ideal_estimate(f_seq ,delta_f ,N ,window_type ,estimate_type)
bias_error = zeros(1,length(f_seq));    %Ԥ����
f_estimate = zeros(1,length(f_seq));    %Ԥ����
for n = 1:length(f_seq)
    t = (0:N-1)/(N*delta_f);    % ����ʱ������
    A = 1;                      % �źŷ���
    x = A*cos(2*pi*f_seq(n)*t);          % �źŲ�������
   %% �Ӵ�
   switch(window_type)
       case('none')
            windows = ones(1,length(x));   %���Ӵ���������ȫ1����
       case('hanning')
            windows = hann(length(x))';   %������
       case('hamming')
            windows = hamming(length(x))';    %������
       case('blackman')        
            windows = blackman(length(x))';   %����������
       case('blackmanharris')
            windows = blackmanharris(length(x))';  %��������-����˹��
   end
   signal_window = x.*windows;        %�źų��ϴ�����
   %% ����
    switch(estimate_type)
        case('direct')
            f_estimate(n) = direct_estimate(signal_window ,delta_f ,N); %ֱ�ӹ���
        case('quinn&AM')
            f_estimate(n) = quinn_AM(signal_window,N,N*delta_f);    %quinn-A&M�㷨����
        case('jacobsen')
            f_estimate(n) = Jacobsen_Interpolation(N*delta_f,signal_window,window_type);    %jacobsen��ֵ����
        case('segment FFT')
            f_estimate(n) = Segment_FFT(N*delta_f,signal_window,pi/10);      %�ֶ�FFT�㷨����
    end
    bias_error(n) = f_estimate(n) - f_seq(n);      % ����bias
end
