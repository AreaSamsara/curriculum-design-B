function [mean_f_estimate ,RMSerror] = fre_estimate(f ,delta_f ,SNR_dB ,N ,N_Iteration ,noise_type ,window_type ,estimate_type)
%����������[����Ƶ�ʾ�ֵ��RMS���] = noise_estimate(�ź�Ƶ�� ,Ƶ�ʷֱ��� ,dB����� ,FFT���� ,�������� ,���������� ,����������)
    for i = 1:length(f)
        %% ��ʼ��
        t = (0:N-1)/(N*delta_f);    %����ʱ������
        signal = cos(2*pi*f(i)*t);  %�����ź�

        %% Ԥ����
        error_iteration = zeros(1,N_Iteration);
        f_iteration = zeros(1,N_Iteration);
        
        for x = 1:N_Iteration
            %% ������
            Eb = 0.5;       %���Ҳ�ƽ������=0.5*(����^2),���з���Ϊ1
            Eb_N0 = 10.^(SNR_dB/10);    %��dB�����ת��Ϊ���������
            N0 = Eb/Eb_N0;      %������������
            [signal_noise] = noise_add(signal,N0,noise_type);   %���źŵ�������
            
            %% �Ӵ�
            switch(window_type)
            case('none')
                windows = ones(1,length(signal));       %���Ӵ���������ȫ1����
            case('hanning')
                windows = hann(length(signal))';       %������
            case('hamming')
                windows = hamming(length(signal))';    %������
            case('blackman')        
                windows = blackman(length(signal))';   %����������
            case('blackmanharris')
                windows = blackmanharris(length(signal))';  %��������-����˹��
            end
            signal_window = signal_noise.*windows;      %�źų��ϴ�����

            %% ����
            switch(estimate_type)
            case('direct')
                f_estimate = direct_estimate(signal_window ,delta_f ,N);    %ֱ�ӹ���

            case('quinn&AM')
                f_estimate = quinn_AM(signal_window,N,N*delta_f);   %quinn-A&M�㷨����

            case('jacobsen')
                f_estimate = Jacobsen_Interpolation(N*delta_f,signal_window,window_type);   %jacobsen��ֵ����

            case('segment FFT')
                f_estimate = Segment_FFT(N*delta_f,signal_window,pi/3); %�ֶ�FFT�㷨����
                
            end
            
        
            error_iteration(1,x) = abs(f(i) - f_estimate);  %��bias
            f_iteration(1,x) = f_estimate;      %��¼���ƽ��
            
            %% ������ʾ
            if(rem(x,10) == 0)
               clc;
               disp('RMSerror_frequency'); 
               disp([estimate_type]);       %��ʾ����������
               disp(['����:',num2str((x+(i-1)*N_Iteration)/...
               (N_Iteration*length(f))*100),'%']);      %��ʾ���Ȱٷֱ�
            end
        end

    %% �� RMSƵ����� �� ƽ������Ƶ��
     RMSerror(i) = sqrt(mean(error_iteration.^2));       %��RMS���
     mean_f_estimate(i) = mean(f_iteration);            %��ƽ�����ƽ��
    end
