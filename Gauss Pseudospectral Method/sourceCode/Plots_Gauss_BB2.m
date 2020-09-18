%--------------------------------------------------------------------------
% Plots_Gauss_BB2.m
% Plotting function, Gauss pseudospectral method
% For Bang-Bang problem, using two-phase approach
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Underlying code: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
function Plots_Gauss_BB2(yA,vA,uA,yB,vB,uB,Lam,p,method)
    %----------------------------------------------------------------------
    fontlabel = 20; % x,y label font size
    fontlegend = 12; % x,y legend font size
    fonttick = 12; % x,y rick font size
    wcolor = [1 1 1]; % white color
    bcolor = [0 0 0]; % black color
    mycolor = lines(8);
    
    set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
    set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
    set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
    
    TA = linspace(p.t0,p.ts,1000);
    TB = linspace(p.ts,p.tf,1000);
    T = [TA,TB(2:end)];
    longtA = cat(1, p.t0, p.tA);
    longtB = cat(1, p.ts, p.tB);
    longt = [longtA; longtB];
    longertA = cat(1, longtA, p.ts);
    longertB = cat(1, longtB, p.tf);
    longert = [longt; p.tf];
    
    Lam1A = Lam.Lam1A; Lam2A = Lam.Lam2A;
    Lam1B = Lam.Lam1B; Lam2B = Lam.Lam2B;
    muA = Lam.muA; muB = Lam.muB;
    
    if strcmp(method,'Pseudospectral')
        ylA = LagrangeInter(longtA',yA',TA); 
        vlA = LagrangeInter(longtA',vA',TA);
        ulA = LagrangeInter(p.tA',uA',TA);
        Lam1lA = LagrangeInter(longertA',Lam1A',TA);
        Lam2lA = LagrangeInter(longertA',Lam2A',TA);
        mulA = LagrangeInter(p.tA',muA',TA);
        ylB = LagrangeInter(longtB',yB',TB); 
        vlB = LagrangeInter(longtB',vB',TB);
        ulB = LagrangeInter(p.tB',uB',TB);
        Lam1lB = LagrangeInter(longertB',Lam1B',TB);
        Lam2lB = LagrangeInter(longertB',Lam2B',TB);
        mulB = LagrangeInter(p.tB',muB',TB);
        yl = [ylA, ylB(2:end)];
        vl = [vlA, vlB(2:end)];
        ul = [ulA, ulB(2:end)];
        Lam1l = [Lam1lA, Lam1lB(2:end)];
        Lam2l = [Lam2lA, Lam2lB(2:end)];
        mul = [mulA, mulB(2:end)];
    end
    
    y = [yA; yB];
    v = [vA; vB];
    u = [uA; uB];
    Lam1 = [Lam1A; Lam1B(2:end)];
    Lam2 = [Lam2A; Lam2B(2:end)];
        
    %----------------------------------------------------------------------
    % plot states
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
 
    [Xopt1, Xopt2] = BangBang_Solution_States(T);
    plot(T,Xopt1,'-','linewidth',2,'color',mycolor(6,:)); hold on
    plot(T,Xopt2,'-','linewidth',2,'color',mycolor(2,:)); hold on
        
    plot(longt,y,'.','markersize',14,'color',mycolor(1,:)); hold on
    plot(longt,v,'.','markersize',14,'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,yl,'-','color',mycolor(1,:)); hold on
        plot(T,vl,'-','color',mycolor(7,:)); hold on
    end
    
    mytitle = ['Gauss ',method,' Method - States']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'states'; % y label with latex
    mylegend = {'True Solution - State 1','True Solution - State 2',...
        [method, ' - State 1'],[method, ' - State 2']}; % legend with latex
    xlabel(myxlabel) % create x label
    ylabel(myylabel) % create y label
    ha = gca; % get current axis handle
    try
        ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
        ha.XAxis.FontSize = fonttick; % change x tick font size
        ha.YAxis.FontSize = fonttick; % change y tick font size
        ha.XAxis.Label.FontSize = fontlabel; % change x label font size
        ha.YAxis.Label.FontSize = fontlabel; % change y label font size
        ht = title(mytitle);
        ht.FontSize = fontlabel;
        hl = legend(mylegend,'location','Best'); % create legend
        hl.FontSize = fontlegend; % change legend font size
        hl.EdgeColor = bcolor; % change the legend border to black (not a dark grey)
    catch
        disp('plot formatting failed (try using a version that supports HG2)')
    end
    
    
    %----------------------------------------------------------------------
    % plot costates
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
 
    [Xopt1, Xopt2] = BangBang_Solution_Costates(T);
    plot(T,Xopt1,'-','linewidth',2,'color',mycolor(6,:)); hold on
    plot(T,Xopt2,'-','linewidth',2,'color',mycolor(2,:)); hold on
        
    plot(longert,Lam1,'.','markersize',14,'color',mycolor(1,:)); hold on
    plot(longert,Lam2,'.','markersize',14,'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,Lam1l,'-','color',mycolor(1,:)); hold on
        plot(T,Lam2l,'-','color',mycolor(7,:)); hold on
    end
    
    mytitle = ['Gauss ',method,' Method - Costates']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'costates'; % y label with latex
    mylegend = {'True Solution - Costate 1','True Solution - Costate 2',...
        [method, ' - Costate 1'],[method, ' - Costate 2']}; % legend with latex
    xlabel(myxlabel) % create x label
    ylabel(myylabel) % create y label
    ha = gca; % get current axis handle
    try
        ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
        ha.XAxis.FontSize = fonttick; % change x tick font size
        ha.YAxis.FontSize = fonttick; % change y tick font size
        ha.XAxis.Label.FontSize = fontlabel; % change x label font size
        ha.YAxis.Label.FontSize = fontlabel; % change y label font size
        ht = title(mytitle);
        ht.FontSize = fontlabel;
        hl = legend(mylegend,'location','Best'); % create legend
        hl.FontSize = fontlegend; % change legend font size
        hl.EdgeColor = bcolor; % change the legend border to black (not a dark grey)
    catch
        disp('plot formatting failed (try using a version that supports HG2)')
    end
    
    
    %----------------------------------------------------------------------
    % plot control
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
    
    Uopt = BangBang_Solution_Control(T);
    plot(T,Uopt,'linewidth',1,'color',mycolor(6,:)); hold on
    
    plot(p.t,u,'.','markersize',14,'color',mycolor(1,:)); hold on

    if strcmp(method,'Pseudospectral')
        plot(T,ul,'-','color',mycolor(1,:)); hold on
    end
    
	mytitle = ['Gauss ',method,' Method - Control']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'control'; % y label with latex
    mylegend = {'True Solution - Control',[method, ' - Control']}; % legend with latex
    xlabel(myxlabel) % create x label
    ylabel(myylabel) % create y label
    try
        ha = gca; % get current axis handle
        ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
        ha.XAxis.FontSize = fonttick; % change x tick font size
        ha.YAxis.FontSize = fonttick; % change y tick font size
        ha.XAxis.Label.FontSize = fontlabel; % change x label font size
        ha.YAxis.Label.FontSize = fontlabel; % change y label font size
        ht = title(mytitle);
        ht.FontSize = fontlabel;
        hl = legend(mylegend,'location','Best'); % create legend
        hl.FontSize = fontlegend; % change legend font size
        hl.EdgeColor = bcolor; % change the legend border to black (not a dark grey)
    catch
        disp('plot formatting failed (try using a version that supports HG2)')
    end
    
    %----------------------------------------------------------------------
    % Calculate and plot Hamiltonian
    Ham = 1 + Lam1l.*vl + Lam2l.*ul - mul.*(abs(ul)-1);
    
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color

    if strcmp(method,'Pseudospectral')
        plot(T,Ham,'-','color',mycolor(1,:),'linewidth',1); hold on
    end
    
	mytitle = ['Gauss ',method,' - Hamiltonian']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'Hamiltonian'; % y label with latex
    ylim([-.00001, .00001])
    xlabel(myxlabel) % create x label
    ylabel(myylabel) % create y label
    try
        ha = gca; % get current axis handle
        ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
        ha.XAxis.FontSize = fonttick; % change x tick font size
        ha.YAxis.FontSize = fonttick; % change y tick font size
        ha.XAxis.Label.FontSize = fontlabel; % change x label font size
        ha.YAxis.Label.FontSize = fontlabel; % change y label font size
        ht = title(mytitle);
        ht.FontSize = fontlabel;
    catch
        disp('plot formatting failed (try using a version that supports HG2)')
    end
end