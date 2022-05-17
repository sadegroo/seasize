function plotpat_multi_hiponly(pat,H_all, idvec)
    
    no = length(idvec);
    noofplots = 4;
    [~,pnts] = size(H_all(1).a);
    figure('Name','patfig')
    
    subplot(noofplots,1,1);
    hold on
    for i = 1:no
        plot(1:pnts,H_all(idvec(i)).a);
    end
    xlim([1 101])
    hold off
    grid on
    title("Hip flexion-extension kinetics")
    ylabel({'Angle', '(°)'})
    legend(cellstr(pat(idvec)), 'Location', 'best');
    
    subplot(noofplots,1,2);
    hold on
    for i = 1:no
        plot(1:pnts,H_all(idvec(i)).ad);
    end
    hold off
    xlim([1 101])
    grid on
    ylabel({'Velocity', '(rad/s)'})
    
    subplot(noofplots,1,3);
    hold on
    for i = 1:no
        plot(1:pnts,H_all(idvec(i)).add);
    end
    hold off
    xlim([1 101])
    grid on
    ylabel({'Acceleration' '(rad/s^2)'})
    
    subplot(noofplots,1,4);
    hold on
    for i = 1:no
        plot(1:pnts,-1*H_all(idvec(i)).m_abs);
    end
    hold off
    xlim([1 101])
    grid on
    ylabel({'Joint moment' '(Nm)'})
    xlabel('Gait %')
    
%     subplot(noofplots,1,5);
%     hold on
%     for i = 1:no
%         plot(1:pnts,H_all(idvec(i)).p_abs);
%     end
%     hold off
%     xlim([1 101])
%     grid on
%     ylabel({'Power' '(W)'})
%     xlabel('Gait %')

    sgtitle("Hip flexion-extension kinetics")
    
    shg
    
    %savefig(gcf,[num2str(id) '_' pat{id} '.fig'])
end