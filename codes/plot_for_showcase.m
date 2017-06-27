function plot_for_showcase(protein,stab_run,tstart,tend)

    timevec = protein.stability_runs(stab_run).time_vector;
    numCom = protein.stability_runs(stab_run).num_communities(timevec >= tstart & timevec <= tend); 
    stabval = protein.stability_runs(stab_run).stability(timevec >= tstart & timevec <= tend); 
    time_vector = protein.stability_runs(stab_run).time_vector(timevec >= tstart & timevec <= tend);
    
    set(0,'DefaultAxesFontSize',5);    
    set(0,'defaultlinelinewidth',1);
    
    
    % plot the stability and number of communities graph
    h = figure();
    
    set(h,'visible','off');
    ax = plotyy(time_vector,numCom,time_vector,stabval);
    
    set(ax(1),'YScale','log');
    set(ax(2),'YScale','log');
    set(ax(1),'YTick',[10^0 10^1 10^2 10^3 10^4]);
    set(ax(2),'YTick',[0 0.0001 0.001 0.01 0.1 1]);
    
    set(get(ax(1),'Ylabel'),'String','Number of communities');
    set(get(ax(2),'Ylabel'),'String','Stability');
    
    set(ax(1),'XLim', [time_vector(1) time_vector(end)], ...
        'YLim', [1 10^4], 'XScale','log');
    set(ax(2),'XLim', [time_vector(1) time_vector(end)],'XScale','log');
    
    xlabel('Markov time');
    ylabel('Number of communities');
    
    hold on
    
    for i = 1:length(time_vector)
        % add the vertical line
        lh = line([time_vector(i) time_vector(i)], [1 10^4], 'Color',[0 0.5 0.5]);

        % save the graph
        % 150 pixels/inch by default
        set(h,'PaperUnits','inches','PaperPosition',[0 0 4 4]); 
        saveas(h,[num2str(i-1) '.png'],'png');
        
        delete(lh);
        disp(['Step ' num2str(i) ' done']);
    end
    
end
