function plotpat(pat,H,K,A,id)

    [~,pnts] = size(H.a);
    figure('Name','patfig')
    
    subplot(4,3,1);
    plot(1:pnts,H.a);
    xlim([1 101])
    grid on
    ylabel({'Angle', '(°)'})
    title('Hip')
    
    subplot(4,3,2);
    plot(1:pnts,K.a);
    xlim([1 101])
    grid on
    title('Knee')
    
    subplot(4,3,3);
    plot(1:pnts,A.a);
    xlim([1 101])
    grid on
    title('Ankle')
    
    subplot(4,3,4);
    plot(1:pnts,H.ad);
    xlim([1 101])
    grid on
    ylabel({'Velocity', '(rad/s)'})
    
    subplot(4,3,5);
    plot(1:pnts,K.ad);
    xlim([1 101])
    grid on
    
    subplot(4,3,6);
    plot(1:pnts,A.ad);
    xlim([1 101])
    grid on
    
    subplot(4,3,7);
    plot(1:pnts,H.add);
    xlim([1 101])
    grid on
    ylabel({'Acceleration' '(rad/s^2)'})
    
    subplot(4,3,8);
    plot(1:pnts,K.add);
    xlim([1 101])
    grid on
    
    subplot(4,3,9);
    plot(1:pnts,A.add);
    xlim([1 101])
    grid on
    
    subplot(4,3,10);
    plot(1:pnts,H.m_abs);
    xlim([1 101])
    grid on
    ylabel({'Load moment' '(Nm)'})
    xlabel('Gait %')
    
    subplot(4,3,11);
    plot(1:pnts,K.m_abs);
    xlim([1 101])
    xlabel('Gait %')
    grid on
    
    subplot(4,3,12);
    plot(1:pnts,A.m_abs);
    xlim([1 101])
    xlabel('Gait %')
    grid on

    sgtitle(pat{id})
    
    shg
    
    %savefig(gcf,[num2str(id) '_' pat{id} '.fig'])
end