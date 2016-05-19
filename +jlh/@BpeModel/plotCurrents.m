function obj = plotCurrents(obj)
    % plot anodic and cathodic currents
    %

    %V = (0.05:0.01:0.25);
    V = ((obj.E_m-0.1):0.01:(obj.E_m+0.1));
    % logI0IronOxidation = log(i0IronOxidation(V));
    % logI0OxygenReduction = log(-i0OxygenReduction(V));

    scrsz = get(0,'ScreenSize');
    f = figure('Position',[1 1 scrsz(3)/3 scrsz(4)],'visible',obj.plotsVisible);
    %hf = figure(1);
     j = 4;
     k = 1;
    h1 = subplot(j,k,1);
    h2 = subplot(j,k,2);
    h3 = subplot(j,k,3);
    h4 = subplot(j,k,4);
    hold(h1,'all');
    %set(gca,'LineStyleOrder', 'r-|g:|b--|k-.')
    lineStyles = {'r--','r:','r-','g--','g:','g-','b--','b:','b-','k--','k:','k-','y--','y:','y-'};
    %nLineStyles = size(lineStyles,2);
    %lg = cell(nReactions*3,1);
    lg = {};
    pi = 1; % plot index
    for j=1:obj.nReactions
        li1 = mod((j-1)*3,size(lineStyles,2))+1;
        li2 = mod((j-1)*3+1,size(lineStyles,2))+1;
        li3 = mod((j-1)*3+2,size(lineStyles,2))+1;
    %     plot(h1,    V,i_cathodic{j}(V),lineStyles{li1}, ...
    %                 V,i_anodic{j}(V), lineStyles{li2}, ...
    %                 V,i{j}(V), lineStyles{li3});
    if obj.i_cathodic{j}(V) ~= 0
        plot(h1, V,obj.i_cathodic{j}(V),lineStyles{li1} );
        lg{pi} = strcat(obj.reactions{j}.name,' cathodic current');
        pi = pi+1;
    end
    if obj.i_anodic{j}(V) ~= 0
        plot(h1, V,obj.i_anodic{j}(V), lineStyles{li2});
        lg{pi} = strcat(obj.reactions{j}.name,' anodic current');
        pi = pi+1;
    end
        plot(h1, V,obj.i{j}(V), lineStyles{li3});
        lg{pi} = strcat(obj.reactions{j}.name,' current');
        pi = pi+1;
    end
    legend(h1,lg,'Location','EastOutside','FontSize',6);
    title(h1,'reaction currents');
    xlabel(h1,'E / V');
    ylabel(h1,'i / A m^{-2}');
    hold(h1,'off');

    hold(h2,'all');
    lineStyles = {'r-','r:','g-','g:','b-','b:','k-','k:','y-','y:'};
    lg = {};
    pi = 1;
    for j=1:obj.nReactions
        li1 = mod((j-1)*2,size(lineStyles,2))+1;
        li2 = mod((j-1)*2+1,size(lineStyles,2))+1;
        if obj.i_cathodic{j}(V) ~= 0
            plot(h2,    V,log(obj.i_cathodic{j}(V)),lineStyles{li1} );
            lg{pi} = strcat(obj.reactions{j}.name,' cathodic current');
            pi = pi+1;
        end
        if obj.i_anodic{j}(V) ~= 0
            plot(h2, V,log(obj.i_anodic{j}(V)),lineStyles{li2});
            lg{pi} = strcat(obj.reactions{j}.name,' anodic current');
            pi = pi+1;
        end
    end
    legend(h2,lg,'Location','EastOutside','FontSize',6);
    title(h2,'log currents');
    xlabel(h2,'E / V');
    ylabel(h2,'log ( i / A m^{-2})');
    hold(h2,'off');

    lg = {'cathodic', 'anodic', 'total'};
    plot(h3,V,obj.i_cathodic_tot(V),V,obj.i_anodic_tot(V),V,obj.i_tot(V));
    legend(h3,lg,'Location','EastOutside','FontSize',6);
    title(h3,'overall currents');
    xlabel(h3,'E / V');
    ylabel(h3,'i / A m^{-2}');

    %plot(h3,V,i_cathodic_tot_log(V),V,-i_anodic_tot_log(V),V,i_tot_log(V));
    lg = {'cathodic', 'anodic'};
    plot(h4,V,log(obj.i_cathodic_tot(V)),V,log(obj.i_anodic_tot(V)));
    legend(h4,lg,'Location','EastOutside','FontSize',6);
    title(h4,'overall currents');
    xlabel(h4,'E / V');
    ylabel(h4,'log (i / A m^{-2})');

    obj.savePlot(f,'matlabCurrentPlots.png',4,2);
end