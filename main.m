% 交互界面demo【如果退不出去了按ctrl+C】
% 信号频率——situation.f_seq
% 干扰频率——situation.Interfere
% 噪声类型——situation.noise_type
% 迭代次数——situation.Iteration
% 信噪比/信干比——situation.SNR_dB
% 选择窗函数，直接修改windows_type
% 选择估计方法，直接修改estimate_method

addpath(genpath(pwd));
clear;close all;
%% 交互界面
stop = false;

while(stop == false)        
    % interactiveUI();
    % break;

    % sel = menu('演示界面菜单',...
    %     '无噪声：估计偏差-信号频率',...
    %     '有噪声：高信噪比下，RMS误差-信号频率',...
    %     '有噪声：低信噪比下，RMS误差-信号频率',...
    %     '有噪声：RMS误差-信噪比',...
    %     '有噪声：不同噪声类型，RMS误差-信噪比',...
    %     '有噪声：不同窗函数，RMS误差-信噪比',...
    %     '有干扰：高信干比下，估计偏差-干扰相对位置',...
    %     '有干扰：低信干比下，估计偏差-干扰相对位置',...
    %     '有干扰：估计偏差-信干比',...
    %     '自定义模式',...
    %     'Exit',...
    %     '窗口');

    % 替换原始的menu函数为更美观的listdlg
    options = {
        '无噪声：估计偏差-信号频率',...
        '有噪声：高信噪比下，RMS误差-信号频率',...
        '有噪声：低信噪比下，RMS误差-信号频率',...
        '有噪声：RMS误差-信噪比',...
        '有噪声：不同噪声类型，RMS误差-信噪比',...
        '有噪声：不同窗函数，RMS误差-信噪比',...
        '有干扰：高信干比下，估计偏差-干扰相对位置',...
        '有干扰：低信干比下，估计偏差-干扰相对位置',...
        '有干扰：估计偏差-信干比',...
        '自定义模式',...
        'Exit',...
        '窗口'};
    
    [sel, ok] = listdlg(...
        'PromptString', '请选择演示模式 (可单选)',...
        'ListString', options,...
        'SelectionMode', 'single',...
        'ListSize', [300 200],...  % 调整弹窗大小
        'Name', '演示界面菜单');
    
    if ~ok  % 如果用户取消选择
        flag = 0;
        break;
    end

    stop = runEstimate(sel);
    if stop
        break;
    end

    % 添加等待环节
    h = uicontrol('Style', 'pushbutton', 'String', '返回菜单',...
                 'Position', [10 10 100 30],...
                 'Callback', 'uiresume(gcbf)');
    uiwait(gcf);  % 暂停直到点击按钮
    delete(h);    % 删除临时按钮

end

function stop = runEstimate(sel)
    %% 初始化参数
    fs = 2560;  % 采样频率
    N = 256;    % 采样点数
    
    stop = false;

    switch(sel)
        case(1)
            %% 无噪声：估计偏差-信号频率
            situation.name = 'ideal';
            situation.f_seq = 210:0.5:220;  %设置仿真频率，需要时可以修改
            situation.SNR_dB = [];
            situation.Iteration = [];
            situation.noise_type = [];
            situation.Interfere = [];
            % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};    % 选择窗函数，可修改
            estimate_method = {'direct','quinn&AM','jacobsen','segment FFT'};   %选择估计算法，可修改
            % estimate_method = {'jacobsen','segment FFT'};
            % estimate_method = {'jacobsen'};
            simlus(fs,N,situation,windows_type,estimate_method);
            case(2)
       %% 有噪声：高信噪比下，RMS误差-信号频率
            situation.name = 'noise';
            situation.f_seq = 210:0.5:220;  % 设置仿真频率，可修改
            situation.SNR_dB = 10;          % 设置信噪比，可改，但只能是一个值
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};      % 噪声类型，可改，但只能为一个
            windows_type = {'none'};                % 窗函数类型，可改，但只能为一个
            estimate_method = {'direct','quinn&AM','jacobsen','segment FFT'};   % 估计算法，可改，但必须为多个
            simlus(fs,N,situation,windows_type,estimate_method);
        case(3)
        %% 有噪声：低信噪比下，RMS误差-信号频率
            situation.name = 'noise';
            situation.f_seq = 210:0.5:220;  % 设置仿真频率，可修改
            situation.SNR_dB = 0;           % 设置信噪比，可改，但只能是一个值
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};       % 噪声类型，可改，但只能为一个
            windows_type = {'none'};                 % 窗函数类型，可改，但只能为一个
            estimate_method = {'direct','quinn&AM','jacobsen','segment FFT'};   % 估计算法，可改，但必须为多个
            simlus(fs,N,situation,windows_type,estimate_method);
        case(4)
       %% 有噪声：RMS误差-信噪比
            situation.name = 'noise';   
            situation.f_seq = 215;          % 设置仿真频率，可改,但只能是一个值
            situation.SNR_dB = 0:0.5:10;    % 设置信噪比，可改
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};  % 噪声类型，可改，但只能为一个
            windows_type = {'none'};            % 窗函数类型，可改，但只能为一个
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};    % 估计算法，可改，但必须为多个
            simlus(fs,N,situation,windows_type,estimate_method);
        case(5)
       %% 有噪声：不同噪声类型，RMS误差-信噪比
            situation.name = 'noise';
            situation.f_seq = 215;          % 设置仿真频率，可改,但只能是一个值
            situation.SNR_dB = 0:0.5:10;    % 设置信噪比，可改
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal','Poisson','Chisquare'};    % 噪声类型，可改，但必须为多个
            windows_type = {'none'};                    % 窗函数类型，可改，但只能为一个
            estimate_method = {'segment FFT'};          % 估计算法，可改，但只能为一个
            simlus(fs,N,situation,windows_type,estimate_method);
        case(6)
       %% 有噪声：不同窗函数，RMS误差-信噪比
            situation.name = 'noise';
            situation.f_seq = 220;          % 设置仿真频率，可改,但只能是一个值
            situation.SNR_dB = 0:0.5:10;    % 设置信噪比，可改
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};      % 噪声类型，可改，但只能为一个
            windows_type = {'none','hanning','hamming','blackman','blackmanharris'};    % 窗函数类型，可改，但必须为多个
            %estimate_method = {'quinn&AM'};         % 估计算法，可改，但只能为一个
            %estimate_method = {'jacobsen'};
            estimate_method = {'segment FFT'};
            simlus(fs,N,situation,windows_type,estimate_method);
        case(7)
       %% 有干扰：高信干比下，估计偏差-干扰相对位置
            situation.name = 'Interfere';
            situation.f_seq = 215;      % 设置仿真频率，可改,但只能是一个值
            situation.SNR_dB = 10;      % 设置信干比，可改，但只能是一个值
            situation.Iteration = [];   
            situation.noise_type = [];
            situation.Interfere = 210:0.5:220;  %设置单频干扰频率
    % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};            % 窗函数类型，可改
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};    % 估计算法，可改
            simlus(fs,N,situation,windows_type,estimate_method);
        case(8)
       %% 有干扰：低信干比下，估计偏差-干扰相对位置
            situation.name = 'Interfere';
            situation.f_seq = 215;      % 设置仿真频率，可改,但只能是一个值
            situation.SNR_dB = 0;       % 设置信干比，可改，但只能是一个值
            situation.Iteration = [];
            situation.noise_type = [];
            situation.Interfere = 210:0.5:220;      %设置单频干扰频率
    % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};            % 窗函数类型，可改
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};    % 估计算法，可改
            simlus(fs,N,situation,windows_type,estimate_method);
        case(9)
       %% 有干扰：估计偏差-信干比
            situation.name = 'Interfere';
            situation.f_seq = 215;      % 设置仿真频率，可改,但只能是一个值
            situation.SNR_dB = 0:0.5:10;    % 设置信干比，可改
            situation.Iteration = [];
            situation.noise_type = [];
            situation.Interfere = 220;      %设置单频干扰频率，可改,但只能是一个值
    % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};            % 窗函数类型，可改
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};        % 估计算法，可改
            simlus(fs,N,situation,windows_type,estimate_method);
        case(10)
       %% 自定义模式
            disp("理想条件'ideal'/有噪声条件'noise'/有干扰条件'Interfere'")
            situation.name = input('仿真条件：');
            if(strcmp(situation.name,'ideal'))
           %% 理想情况
                situation.f_seq = 210:0.5:220;  % 设置仿真频率，可改
                situation.SNR_dB = [];
                situation.Iteration = [];
                situation.noise_type = [];
                situation.Interfere = [];
                disp("不加窗(矩形窗)'none'/汉宁窗'hanning'/汉明窗'hamming'/布莱克曼窗'blackman'/布莱克曼-哈里斯窗'blackmanharris'")
                windows_type = input('选择窗函数(用元胞表示)：');
                disp("直接法'direct'/A&M迭代联合Quinn法'quinn&AM'/Jacobsen改进插值法'jacobsen'/分段FFT法'segment FFT'")
                estimate_method = input('选择估计算法(用元胞表示)：');
            elseif(strcmp(situation.name,'Interfere'))
           %% 有干扰情况
                disp("信干比'SIR'/相对位置'df'")
                mod_sel = input('选择横坐标：');
                if(strcmp(mod_sel,'SIR'))
              %% 考虑估计偏差-SIR
                    situation.f_seq = 215;          % 设置仿真频率，可改,但只能是一个值
                    situation.SNR_dB = 0:0.5:10;    % 设置信干比，可改
                    situation.Iteration = [];
                    situation.noise_type = [];
                    situation.Interfere = 220;      % 设置单频干扰频率，可改，但只能是一个值
                    disp("不加窗(矩形窗)'none'/汉宁窗'hanning'/汉明窗'hamming'/布莱克曼窗'blackman'/布莱克曼-哈里斯窗'blackmanharris'")
                    windows_type = input('选择窗函数(用元胞表示)：');
                    disp("直接法'direct'/A&M迭代联合Quinn法'quinn&AM'/Jacobsen改进插值法'jacobsen'/分段FFT法'segment FFT'")
                    estimate_method = input('选择估计算法(用元胞表示)：');
                else
              %% 考虑估计偏差-相对位置
                    situation.name = 'Interfere';   
                    situation.f_seq = 215;      % 设置仿真频率，可改,但只能是一个值
                    situation.SNR_dB = input('信干比为：');
                    situation.Iteration = [];
                    situation.noise_type = [];
                    situation.Interfere = 210:0.5:220;  % 设置单频干扰频率，可改
                    disp("不加窗(矩形窗)'none'/汉宁窗'hanning'/汉明窗'hamming'/布莱克曼窗'blackman'/布莱克曼-哈里斯窗'blackmanharris'")
                    windows_type = input('选择窗函数(用元胞表示)：');
                    disp("直接法'direct'/A&M迭代联合Quinn法'quinn&AM'/Jacobsen改进插值法'jacobsen'/分段FFT法'segment FFT'")
                    estimate_method = input('选择估计算法(用元胞表示)：');
                end
            elseif(strcmp(situation.name,'noise'))
           %% 有噪声情况
                disp("信噪比'SNR'/信号频率'f'")
                mod_sel = input('选择横坐标：');
                if(strcmp(mod_sel,'SNR'))
              %% 考虑RMS偏差-信噪比
                    situation.f_seq = 215;       % 设置仿真频率，可改,但只能是一个值
                    situation.SNR_dB = 0:0.5:10;    % 设置信噪比，可改
                    situation.Iteration = 1e3;
                    disp("高斯分布'Normal'/泊松分布'Poisson'/卡方分布'Chisquare'")
                    situation.noise_type = input('选择噪声类型(用元胞表示)：');
                    disp("不加窗(矩形窗)'none'/汉宁窗'hanning'/汉明窗'hamming'/布莱克曼窗'blackman'/布莱克曼-哈里斯窗'blackmanharris'")
                    windows_type = input('选择窗函数(用元胞表示)：');
                    disp("直接法'direct'/A&M迭代联合Quinn法'quinn&AM'/Jacobsen改进插值法'jacobsen'/分段FFT法'segment FFT'")
                    estimate_method = input('选择估计算法(用元胞表示)：');
                else
              %% 考虑RMS偏差-信号频率
                    situation.f_seq = 210:0.5:220;      % 设置仿真频率，可改
                    situation.SNR_dB = input('选择信噪比：(dB)');
                    situation.Iteration = 1e3;
                    disp("高斯分布'Normal'/泊松分布'Poisson'/卡方分布'Chisquare'")
                    situation.noise_type = input('选择噪声类型(用元胞表示)：');
                    disp("不加窗(矩形窗)'none'/汉宁窗'hanning'/汉明窗'hamming'/布莱克曼窗'blackman'/布莱克曼-哈里斯窗'blackmanharris'")
                    windows_type = input('选择窗函数(用元胞表示)：');
                    disp("直接法'direct'/A&M迭代联合Quinn法'quinn&AM'/Jacobsen改进插值法'jacobsen'/分段FFT法'segment FFT'")
                    estimate_method = input('选择估计算法(用元胞表示)：');
                end
            end
            simlus(fs,N,situation,windows_type,estimate_method);
        case(11)
       %% 退出程序
            stop = true;
        case(12)
            window();
    end
end