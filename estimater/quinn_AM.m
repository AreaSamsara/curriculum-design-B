% 函数输入描述：s=信号序列 N=傅里叶变换点数 fs=采样频率
% 函数输出描述：f_quinn=使用基于quinn算法的A&M改进算法的频率估值
function f_quinn = quinn_AM(s,N,fs)
    %% quinn算法
    s_fft = fft(s,N);       % 进行fft运算
    [~,index] = max(abs(s_fft));      % 寻找最大谱峰
    a1 = real(s_fft(index-1)/s_fft(index));    % 利用最大谱峰与左边的谱线计算中间变量
    a2 = real(s_fft(index+1)/s_fft(index));    % 利用最大谱峰与右边的谱线计算中间变量
    deta1 = a1/(1-a1);     % 利用最大谱峰与左边谱线得出的频率误差
    deta2 = a2/(a2-1);     % 利用最大谱峰与右边谱线得出的频率误差
    % 判断使用哪个频率误差
    if deta1>0&&deta2>0    
        deta_1 = deta2;
    else
        deta_1 = deta1;
    end
    
    %% A&M算法
    index2 = index + deta_1;  
    p = [0.5,-0.5];
    ss = [0,0];
    for i = 1:2
        for j = 0:N-1
            ss(i) = ss(i) + s(j+1)*exp(-1j*2*pi*(index2+p(i))*j/N);  % 迭代式
        end
    end
    deta3 = real((ss(1)+ss(2))/(ss(1)-ss(2)))/2;    % 第二次迭代的频率误差

    deta_2 = deta_1 + deta3;      % 总频率误差
    f_quinn = (index + deta_2)*fs/N;   % 对频率进行修正得到频率估值
end
