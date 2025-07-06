function [f_estimate ,error] = int_estimate_f(f ,f_int ,delta_f ,SIR_dB ,N ,window_type ,estimate_type)
%函数描述：[估计频率，误差] = int_estimate(信号频率 ,干扰频率，频率分辨率 ,dB信干比 ,FFT点数 ,噪声类型，窗函数类型 ,估计器类型)
    Eb = 0.5;       %正弦波平均功率=0.5*(幅度^2),其中幅度为1
    SIR = 10^(SIR_dB/10);   %将dB信干比转换为线性信干比
    I = Eb/SIR;             %计算干扰功率
    a_int = sqrt(2*I);      %计算干扰幅度

    t = (0:N-1)/(delta_f*N);    %生成时间序列

    for i = 1:length(f_int)
        %% 加入干扰
        signal_int = cos(2*pi*f*t) + a_int*cos(2*pi*f_int(i)*t);

        %% 加窗
        switch(window_type)
        case('none')
            windows = ones(1,length(signal_int));   %不加窗，即乘上全1矩阵
        case('hanning')
            windows = hann(length(signal_int))';   %汉宁窗
        case('hamming')
            windows = hamming(length(signal_int))';    %汉明窗
        case('blackman')        
            windows = blackman(length(signal_int))';   %布莱克曼窗
        case('blackmanharris')
            windows = blackmanharris(length(signal_int))';  %布莱克曼-哈里斯窗
        end
            
        signal_window = signal_int.*windows;        %信号乘上窗函数


        %% 估计
        switch(estimate_type)
        case('direct')
            f_est = direct_estimate(signal_window ,delta_f ,N); %直接估计
        case('quinn&AM')
            f_est = quinn_AM(signal_window,N,N*delta_f);    %quinn-A&M算法估计
        case('jacobsen')
            f_est = Jacobsen_Interpolation(N*delta_f,signal_window,window_type);    %jacobsen插值估计
        case('segment FFT')
            f_est = Segment_FFT(N*delta_f,signal_window,pi/3);      %分段FFT算法估计
        end 
            
        %% 求 频率误差 和 估计频率
        error(1,i) = abs(f - f_est);        %求bias
        f_estimate(1,i) = f_est;            %记录估计结果
        
       clc;
       disp('error_int_f');
       disp([estimate_type]);               %显示估计器类型
       disp(['进度:',num2str(i/length(f_int)*100),'%']);   %显示进度百分比
    end
