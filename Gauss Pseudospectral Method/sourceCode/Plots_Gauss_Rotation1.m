%--------------------------------------------------------------------------
% Plots_Gauss_Rotation1.m
% Plotting function, Gauss pseudospectral method
% For rotation problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Underlying code: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
function Plots_Gauss_Rotation1(q,w,u,p,method)
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
    
    if strcmp(method,'Pseudospectral')
        q1l = LagrangeInter(longt',q(:,1)',T);
        q2l = LagrangeInter(longt',q(:,2)',T);
        q3l = LagrangeInter(longt',q(:,3)',T);
        q4l = LagrangeInter(longt',q(:,4)',T);
        w1l = LagrangeInter(longt',w(:,1)',T);
        w2l = LagrangeInter(longt',w(:,2)',T);
        w3l = LagrangeInter(longt',w(:,3)',T);
        u1l = LagrangeInter(p.t',u(:,1)',T);
        u2l = LagrangeInter(p.t',u(:,2)',T);
        u3l = LagrangeInter(p.t',u(:,3)',T);
    end
    
    
    %----------------------------------------------------------------------
    % plot quaternion magnitude
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
    
    qMag = sqrt(q1l.^2 + q2l.^2 + q3l.^2 + q4l.^2);
    plot(T,qMag,'-','color',mycolor(7,:)); hold on
    
    mytitle = ['Quaternion Magnitude']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'Magnitude'; % y label with latex
    xlabel(myxlabel) % create x label
    ylabel(myylabel) % create y label
    ylim([0 1.3])
    ha = gca; % get current axis handle
    try
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
        
    %----------------------------------------------------------------------
    % plot quaternion components
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,q(:,1),'.','markersize',14,'color',mycolor(1,:)); hold on
    plot(longt,q(:,2),'.','markersize',14,'color',mycolor(3,:)); hold on
    plot(longt,q(:,3),'.','markersize',14,'color',mycolor(5,:)); hold on
    plot(longt,q(:,4),'.','markersize',14,'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,q1l,'-','color',mycolor(1,:)); hold on
        plot(T,q2l,'-','color',mycolor(3,:)); hold on
        plot(T,q3l,'-','color',mycolor(5,:)); hold on
        plot(T,q4l,'-','color',mycolor(7,:)); hold on
    end
    
    mytitle = ['Quaternion Components']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'Quaternion Components'; % y label with latex
    mylegend = {'q1','q2','q3','q4'}; % legend with latex
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
    % plot angular velocity components
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,w(:,1),'.','markersize',14,'color',mycolor(2,:)); hold on
    plot(longt,w(:,2),'.','markersize',14,'color',mycolor(4,:)); hold on
    plot(longt,w(:,3),'.','markersize',14,'color',mycolor(6,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,w1l,'-','color',mycolor(2,:)); hold on
        plot(T,w2l,'-','color',mycolor(4,:)); hold on
        plot(T,w3l,'-','color',mycolor(6,:)); hold on
    end
    
    mytitle = ['Angular Velocity Components']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'Angular Velocity (rad/sec)'; % y label with latex
    mylegend = {'$\omega_x$','$\omega_y$','$\omega_z$'}; % legend with latex
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
    % plot torque components
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(p.t,u(:,1),'.','markersize',14,'color',mycolor(2,:)); hold on
    plot(p.t,u(:,2),'.','markersize',14,'color',mycolor(4,:)); hold on
    plot(p.t,u(:,3),'.','markersize',14,'color',mycolor(6,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,u1l,'-','color',mycolor(2,:)); hold on
        plot(T,u2l,'-','color',mycolor(4,:)); hold on
        plot(T,u3l,'-','color',mycolor(6,:)); hold on
    end
    
    mytitle = ['Torque Components']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'Torque (N-m)'; % y label with latex
    mylegend = {'$N_x$','$N_y$','$N_z$'}; % legend with latex
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
    
  