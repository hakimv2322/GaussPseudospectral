%--------------------------------------------------------------------------
% Plots_Gauss_LQR.m
% Plotting function, Gauss pseudospectral method
% For LQR problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Underlying code: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
function Plots_Gauss_LQR(y,u,Lam,p,method)
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
    
    T = linspace(p.t0,p.tf,10000);
    longt = cat(1, p.t0, p.t);
    longert = cat(1, longt, p.tf);
    
    if strcmp(method,'Pseudospectral')
        yl = LagrangeInter(longt',y',T); 
        ul = LagrangeInter(p.t',u',T);
        Laml = LagrangeInter(longert',Lam',T);
    end
        
    %----------------------------------------------------------------------
    % plot states
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
 
    Xopt = LQR_Solution_State(T);
    plot(T,Xopt,'-','linewidth',2,'color',mycolor(6,:)); hold on
        
    plot(longt,y,'.','markersize',14,'color',mycolor(1,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,yl,'-','color',mycolor(1,:)); hold on
    end
    
    mytitle = ['Gauss ',method,' Method - States']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'states'; % y label with latex
    mylegend = {'True Solution - State 1',...
        [method, ' - State 1']}; % legend with latex
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
 
    Xopt = LQR_Solution_Costate(T);
    plot(T,Xopt,'-','linewidth',2,'color',mycolor(6,:)); hold on
        
    plot(longert,Lam,'.','markersize',14,'color',mycolor(1,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,Laml,'-','color',mycolor(1,:)); hold on
    end
    
    mytitle = ['Gauss ',method,' Method - Costates']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'costates'; % y label with latex
    mylegend = {'True Solution - Costate 1',...
        [method, ' - Costate 1']}; % legend with latex
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
    
    Uopt = LQR_Solution_Control(T);
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
    g = (yl.^2 + ul.^2)./2;
    f = yl + ul;
    Ham = g + Laml.*f;
    
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color

    if strcmp(method,'Pseudospectral')
        plot(T,Ham,'-','color',mycolor(1,:),'linewidth',1); hold on
    end
    
	mytitle = ['Gauss ',method,' Method - Hamiltonian']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'Hamiltonian'; % y label with latex
    ylim([-.1, .1])
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