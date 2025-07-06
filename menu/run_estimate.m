function stop = run_estimate(sel)   
    stop = false;

    switch(sel)
        case(1)
            % 比较不同窗函数的影响
            compareWindowFunctionsXwz(); 
        case(2) 
            % 比较目标信号落入不同间隔时的性能
            compareTargetPositionPerformanceCyz();
        case(3)
            % 比较不同SNR下的影响
            compareSNRPerformanceXwz();
        case(4)
            % 比较不同干扰大小下的性能
            compareInterferenceAmtipleXwz();
        case(5)
            % 比较不同干扰间隔下的性能
            compareInterferenceIntervalXwz();
        case(6)

        otherwise
            stop = true;
    end
end