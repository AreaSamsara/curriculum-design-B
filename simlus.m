% 实现各种情况下的仿真
% 输入：
% fs：采样频率
% N：采样点数
% situation：决定仿真条件（结构体）
%       .name：无噪声无干扰的理想条件ideal/有噪声条件noise/有干扰条件Interfere
%       .f_seq：需仿真的频率序列
%       .SNR_dB：需仿真的信噪比/信干比；如果是理想条件就设为空
%       .Iteration：蒙特卡洛迭代次数。仅在使用蒙特卡洛法时有值，其余情况为空
%       .noise_type：噪声类型（元胞数组）
%               高斯分布'Normal'/泊松分布'Poisson'/卡方分布'Chisquare'
%       .Interfere：干扰频率
% windows_type：仿真的窗函数类型（元胞数组）
%       不加窗(矩形窗)'none'/汉宁窗'hanning'/汉明窗'hamming'/布莱克曼窗'blackman'/布莱克曼-哈里斯窗'blackmanharris'
% estimate_method：选用的估计算法（元胞数组）
%       直接法'direct'/A&M迭代联合Quinn法'quinn&AM'/Jacobsen改进插值法'jacobsen'/分段FFT法'segment FFT'
function [] = simlus(fs,N,situation,windows_type,estimate_method)
%% 基本参数
delta_f = fs/N;                 %频率分辨率
%% 理想情况
if(strcmp(situation.name,'ideal'))
    f_seq = situation.f_seq;    %待仿真的频率
    N_curve = length(windows_type)*length(estimate_method); %图框中曲线的个数
    legend_str = cell(1,N_curve);    %存储图例
    % 计算各种条件下的偏差，并画图
    figure;
    hold on;    grid on;
    title('无噪声无干扰理想单频信号，估计误差与信号频率的关系')
    xlabel(['信号频率/Hz(频谱分辨率',num2str(delta_f),')'])
    ylabel('估计偏差/Hz')
    for n_wintype = 1:length(windows_type)
        for n_estmethod = 1:length(estimate_method)
            [f_estimate,bias_error] = ideal_estimate(f_seq ,delta_f ,N ,windows_type{n_wintype} ,estimate_method{n_estmethod});
            plot(f_seq,abs(bias_error))
            legend_str{(n_wintype-1)*length(estimate_method)+n_estmethod} = [windows_type{n_wintype},';',estimate_method{n_estmethod}];
        end
    end
    legend(legend_str)
%% 有噪声情况
elseif(strcmp(situation.name,'noise'))
    f_seq = situation.f_seq;    % 待仿真的频率
    SNR_dB = situation.SNR_dB;  % 待仿真的信噪比
    N_Iteration = situation.Iteration;      % 迭代次数
    noise_type = situation.noise_type;      % 噪声类型
    if(length(SNR_dB)==1 && length(windows_type)==1 && length(noise_type)==1)
        % 固定信噪比、窗函数、噪声类型，改变信号频率，考虑四种估计算法
        N_curve = length(estimate_method); %图框中曲线的个数
        legend_str = cell(1,N_curve);    %存储图例
        figure;
        hold on;    grid on
        title(['信噪比为',num2str(SNR_dB),'dB',...
            '窗函数为',windows_type{1},...
            '噪声分类型为',noise_type{1}]);
        xlabel(['信号频率/Hz(频谱分辨率',num2str(delta_f),')']);
        ylabel('RMS误差')
        for n_estmethod = 1:length(estimate_method)
            [f_estimate,RMSerror] = fre_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{1} ,windows_type{1} ,estimate_method{n_estmethod});
            plot(f_seq,RMSerror)
            legend_str{n_estmethod} = [estimate_method{n_estmethod}];
        end
        legend(legend_str)
    elseif(length(f_seq)==1 && length(windows_type)==1 && length(noise_type)==1)
        % 固定信号频率、窗函数、噪声类型，改变信噪比，考虑4种估计算法
        N_curve = length(estimate_method); %图框中曲线的个数
        legend_str = cell(1,N_curve);    %存储图例
        figure;
        hold on;    grid on
        title(['信号频率为',num2str(f_seq),...
            '窗函数为',windows_type{1},...
            '噪声分类型为',noise_type{1}]);
        xlabel('信噪比/dB');
        ylabel('RMS误差')
        for n_estmethod = 1:length(estimate_method)
            [f_estimate,RMSerror] = noise_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{1} ,windows_type{1} ,estimate_method{n_estmethod});
            plot(SNR_dB,RMSerror)
            legend_str{n_estmethod} = [estimate_method{n_estmethod}];
        end
        legend(legend_str);
    elseif(length(f_seq)==1 && length(windows_type)==1 && length(estimate_method)==1)
        % 固定信号频率、窗函数、估计算法，改变信噪比，考虑3种噪声类型
        N_curve = length(noise_type); %图框中曲线的个数
        legend_str = cell(1,N_curve);    %存储图例
        figure;
        hold on;    grid on
        title(['信号频率为',num2str(f_seq),...
            '窗函数为',windows_type{1},...
            '估计算法为',estimate_method{1}]);
        xlabel('信噪比/dB');
        ylabel('RMS误差')
        for n_ntype = 1:length(noise_type)
            [f_estimate,RMSerror] = noise_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{n_ntype} ,windows_type{1} ,estimate_method{1});
            plot(SNR_dB,RMSerror)
            legend_str{n_ntype} = [noise_type{n_ntype}];
        end
        legend(legend_str);
    elseif(length(f_seq)==1 && length(noise_type)==1 && length(estimate_method)==1)
        % 固定信号频率、噪声类型、估计算法，改变信噪比，考虑5种窗函数类型
        N_curve = length(windows_type); %图框中曲线的个数
        legend_str = cell(1,N_curve);    %存储图例
        figure;
        hold on;    grid on
        title(['信号频率为',num2str(f_seq),...
            '噪声分类型为',noise_type{1},...
            '估计算法为',estimate_method{1}]);
        xlabel('信噪比/dB');
        ylabel('RMS误差')
        for n_wintype = 1:length(windows_type)
            [f_estimate,RMSerror] = noise_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{1} ,windows_type{n_wintype} ,estimate_method{1});
            disp(windows_type{n_wintype})
            plot(SNR_dB,RMSerror)
            legend_str{n_wintype} = [windows_type{n_wintype}];
        end
        legend(legend_str);
    end
elseif(strcmp(situation.name,'Interfere'))
    f_seq = situation.f_seq;    % 待仿真的频率
    SIR_dB = situation.SNR_dB;  % 待仿真的信干比
    f_int = situation.Interfere;    % 单频干扰频率
    if(length(SIR_dB)==1 && length(f_seq)==1 )
    % 固定信干比和信号频率，改变干扰的相对位置
        N_curve = length(windows_type)*length(estimate_method); %图框中曲线的个数
        legend_str = cell(1,N_curve);    %存储图例
        figure;
        hold on;    grid on;
        title(['单频干扰信号下，信干比为',num2str(SIR_dB),'dB,信号频率为',num2str(f_seq),'Hz'])
        xlabel(['相对位置/Hz','(频谱分辨率',num2str(delta_f),')'])
        ylabel('估计偏差/Hz')
        for n_wintype = 1:length(windows_type)
            for n_estmethod = 1:length(estimate_method)
                [f_estimate,bias_error] = int_estimate_f(f_seq ,f_int ,delta_f ,SIR_dB ,N ,windows_type{n_wintype} ,estimate_method{n_estmethod});
                plot(f_int-f_seq,bias_error)
                legend_str{(n_wintype-1)*length(estimate_method)+n_estmethod} = [windows_type{n_wintype},';',estimate_method{n_estmethod}];
            end
        end
        legend(legend_str);
    elseif(length(f_int)==1 && length(f_seq)==1)
    % 固定干扰频率和信号频率，改变信干比
        N_curve = length(windows_type)*length(estimate_method); %图框中曲线的个数
        legend_str = cell(1,N_curve);    %存储图例
        figure;
        hold on;    grid on;
        title(['单频干扰信号下，干扰频率为',num2str(f_int),'Hz,信号频率为',num2str(f_seq),'Hz',...
            '(频谱分辨率',num2str(delta_f),')'])
        xlabel('信干比/dB')
        ylabel('估计偏差/Hz')
        for n_wintype = 1:length(windows_type)
            for n_estmethod = 1:length(estimate_method)
                [f_estimate,bias_error] = int_estimate_SIR(f_seq ,f_int ,delta_f ,SIR_dB ,N ,windows_type{n_wintype} ,estimate_method{n_estmethod});
                plot(SIR_dB,bias_error)
                legend_str{(n_wintype-1)*length(estimate_method)+n_estmethod} = [windows_type{n_wintype},';',estimate_method{n_estmethod}];
            end
        end
        legend(legend_str);
        
    end
end
