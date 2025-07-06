function [mean_f_estimate ,RMSerror] = fre_estimate(f ,delta_f ,SNR_dB ,N ,N_Iteration ,noise_type ,window_type ,estimate_type)
%函数描述：[估计频率均值，RMS误差] = noise_estimate(信号频率 ,频率分辨率 ,dB信噪比 ,FFT点数 ,迭代次数 ,窗函数类型 ,估计器类型)
    for i = 1:length(f)
        %% 初始化
        t = (0:N-1)/(N*delta_f);    %生成时间序列
        signal = cos(2*pi*f(i)*t);  %生成信号

        %% 预分配
        error_iteration = zeros(1,N_Iteration);
        f_iteration = zeros(1,N_Iteration);
        
        for x = 1:N_Iteration
            %% 加噪声
            Eb = 0.5;       %正弦波平均功率=0.5*(幅度^2),其中幅度为1
            Eb_N0 = 10.^(SNR_dB/10);    %将dB信噪比转换为线性信噪比
            N0 = Eb/Eb_N0;      %计算噪声功率
            [signal_noise] = noise_add(signal,N0,noise_type);   %给信号叠加噪声
            
            %% 加窗
            switch(window_type)
            case('none')
                windows = ones(1,length(signal));       %不加窗，即乘上全1矩阵
            case('hanning')
                windows = hann(length(signal))';       %汉宁窗
            case('hamming')
                windows = hamming(length(signal))';    %汉明窗
            case('blackman')        
                windows = blackman(length(signal))';   %布莱克曼窗
            case('blackmanharris')
                windows = blackmanharris(length(signal))';  %布莱克曼-哈里斯窗
            end
            signal_window = signal_noise.*windows;      %信号乘上窗函数

            %% 估计
            switch(estimate_type)
            case('direct')
                f_estimate = direct_estimate(signal_window ,delta_f ,N);    %直接估计

            case('quinn&AM')
                f_estimate = quinn_AM(signal_window,N,N*delta_f);   %quinn-A&M算法估计

            case('jacobsen')
                f_estimate = Jacobsen_Interpolation(N*delta_f,signal_window,window_type);   %jacobsen插值估计

            case('segment FFT')
                f_estimate = Segment_FFT(N*delta_f,signal_window,pi/3); %分段FFT算法估计
                
            end
            
        
            error_iteration(1,x) = abs(f(i) - f_estimate);  %求bias
            f_iteration(1,x) = f_estimate;      %记录估计结果
            
            %% 进度显示
            if(rem(x,10) == 0)
               clc;
               disp('RMSerror_frequency'); 
               disp([estimate_type]);       %显示估计器类型
               disp(['进度:',num2str((x+(i-1)*N_Iteration)/...
               (N_Iteration*length(f))*100),'%']);      %显示进度百分比
            end
        end

    %% 求 RMS频率误差 和 平均估计频率
     RMSerror(i) = sqrt(mean(error_iteration.^2));       %求RMS误差
     mean_f_estimate(i) = mean(f_iteration);            %求平均估计结果
    end
