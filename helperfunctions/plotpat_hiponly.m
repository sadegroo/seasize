function plotpat_hiponly(pat,H,id)

    [~,pnts] = size(H.a);
    figure('Name','patfig')
    
    subplot(4,1,1);
    plot(1:pnts,H.a);
    xlim([1 101])
    grid on
    ylabel({'Angle', '(°)'})
    title('Hip')
    
    subplot(4,1,2);
    plot(1:pnts,H.ad);
    xlim([1 101])
    grid on
    ylabel({'Velocity', '(rad/s)'})
    
    subplot(4,1,3);
    plot(1:pnts,H.add);
    xlim([1 101])
    grid on
    ylabel({'Acceleration' '(rad/s^2)'})
    
    subplot(4,1,4);
    plot(1:pnts,H.m_abs);
    xlim([1 101])
    grid on
    ylabel({'Load moment' '(Nm)'})
    xlabel('Gait %')

    sgtitle(pat{id})
    
    shg
    
    %savefig(gcf,[num2str(id) '_' pat{id} '.fig'])
end