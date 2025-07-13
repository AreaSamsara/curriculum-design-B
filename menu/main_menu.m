function selection = main_menu()
    stop = false;

    while(stop == false)  
        options = {
            '比较不同窗函数的影响',...
            '比较目标信号落入不同间隔时的性能',...
            '比较不同SNR下的影响',...
            '比较不同干扰大小下的性能',...
            '比较不同干扰间隔下的性能',...
            '比较噪声/干扰服从不同概率分布函数的性能'};
        
        [sel, ok] = listdlg(...
            'PromptString', '请选择仿真方法',...
            'ListString', options,...
            'SelectionMode', 'single',...
            'ListSize', [400 300],...  % 调整弹窗大小
            'Name', '仿真选项菜单');
        
        % 如果用户取消选择，直接退出菜单
        if ~ok
            break;
        end

        % 根据选择结果运行对应的估计函数
        stop = run_estimate(sel);

        % 遇到了意外选项，直接退出
        if stop
            break;
        end

        menu_pause();
    end
end