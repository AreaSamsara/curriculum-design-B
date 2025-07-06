function [f_estimate ,error] = int_estimate_f(f ,f_int ,delta_f ,SIR_dB ,N ,window_type ,estimate_type)
%����������[����Ƶ�ʣ����] = int_estimate(�ź�Ƶ�� ,����Ƶ�ʣ�Ƶ�ʷֱ��� ,dB�Ÿɱ� ,FFT���� ,�������ͣ����������� ,����������)
    Eb = 0.5;       %���Ҳ�ƽ������=0.5*(����^2),���з���Ϊ1
    SIR = 10^(SIR_dB/10);   %��dB�Ÿɱ�ת��Ϊ�����Ÿɱ�
    I = Eb/SIR;             %������Ź���
    a_int = sqrt(2*I);      %������ŷ���

    t = (0:N-1)/(delta_f*N);    %����ʱ������

    for i = 1:length(f_int)
        %% �������
        signal_int = cos(2*pi*f*t) + a_int*cos(2*pi*f_int(i)*t);

        %% �Ӵ�
        switch(window_type)
        case('none')
            windows = ones(1,length(signal_int));   %���Ӵ���������ȫ1����
        case('hanning')
            windows = hann(length(signal_int))';   %������
        case('hamming')
            windows = hamming(length(signal_int))';    %������
        case('blackman')        
            windows = blackman(length(signal_int))';   %����������
        case('blackmanharris')
            windows = blackmanharris(length(signal_int))';  %��������-����˹��
        end
            
        signal_window = signal_int.*windows;        %�źų��ϴ�����


        %% ����
        switch(estimate_type)
        case('direct')
            f_est = direct_estimate(signal_window ,delta_f ,N); %ֱ�ӹ���
        case('quinn&AM')
            f_est = quinn_AM(signal_window,N,N*delta_f);    %quinn-A&M�㷨����
        case('jacobsen')
            f_est = Jacobsen_Interpolation(N*delta_f,signal_window,window_type);    %jacobsen��ֵ����
        case('segment FFT')
            f_est = Segment_FFT(N*delta_f,signal_window,pi/3);      %�ֶ�FFT�㷨����
        end 
            
        %% �� Ƶ����� �� ����Ƶ��
        error(1,i) = abs(f - f_est);        %��bias
        f_estimate(1,i) = f_est;            %��¼���ƽ��
        
       clc;
       disp('error_int_f');
       disp([estimate_type]);               %��ʾ����������
       disp(['����:',num2str(i/length(f_int)*100),'%']);   %��ʾ���Ȱٷֱ�
    end
