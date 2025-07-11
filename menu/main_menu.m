function selection = main_menu()
    stop = false;

    while(stop == false)  
        % 替换原始的menu函数为更美观的listdlg
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
        
        if ~ok  % 如果用户取消选择
            break;
        end

        stop = run_estimate(sel);
        if stop
            break;
        end

        menu_pause();
    end
end