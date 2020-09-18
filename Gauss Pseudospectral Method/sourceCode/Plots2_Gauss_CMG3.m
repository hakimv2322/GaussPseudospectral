%--------------------------------------------------------------------------
% Plot2s_Gauss_CMG3.m
% Plotting function, Gauss pseudospectral method
% For 3/4 CMG pyramid configuration
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Underlying code: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
function Plots2_Gauss_CMG3(q,w,delta,u,Det,p,Wie,method)
    %----------------------------------------------------------------------
    fontlabel = 24; % x,y label font size
    fontlegend = 18; % x,y legend font size
    fonttick = 18; % x,y rick font size
    wcolor = [1 1 1]; % white color
    bcolor = [0 0 0]; % black color
    mycolor = lines(8);
    
    set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
    set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
    set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
    
    set(0, 'DefaultLineLineWidth', 1.2);
    
    T = linspace(p.t0,p.tf,10000);
    longt = [p.t0; p.t];
%     longert = [longt; p.tf];
    
    if strcmp(method,'Pseudospectral')
        q1l = LagrangeInter(longt',q(:,1)',T);
        q2l = LagrangeInter(longt',q(:,2)',T);
        q3l = LagrangeInter(longt',q(:,3)',T);
        q4l = LagrangeInter(longt',q(:,4)',T);
        w1l = LagrangeInter(longt',w(:,1)',T);
        w2l = LagrangeInter(longt',w(:,2)',T);
        w3l = LagrangeInter(longt',w(:,3)',T);
        d1l = LagrangeInter(longt',delta(:,1)',T);
        d2l = LagrangeInter(longt',delta(:,2)',T);
        d3l = LagrangeInter(longt',delta(:,3)',T);
        u1l = LagrangeInter(p.t',u(:,1)',T);
        u2l = LagrangeInter(p.t',u(:,2)',T);
        u3l = LagrangeInter(p.t',u(:,3)',T);
        Detl = LagrangeInter(longt',Det',T);
    end
    
    
    %----------------------------------------------------------------------
    % Plot Q1
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,q(:,1),'.','markersize',14,'color',mycolor(1,:)); hold on
    
    plot(Wie.tArray,Wie.qArray(1,:),'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,q1l,'-','color',mycolor(1,:)); hold on
    end
    
    mytitle = ['$q_1$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = ''; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot Q2
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,q(:,2),'.','markersize',14,'color',mycolor(1,:)); hold on
    
    plot(Wie.tArray,Wie.qArray(2,:),'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,q2l,'-','color',mycolor(1,:)); hold on
    end
    
    mytitle = ['$q_2$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = ''; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot Q3
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,q(:,3),'.','markersize',14,'color',mycolor(1,:)); hold on
    
    plot(Wie.tArray,Wie.qArray(3,:),'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,q3l,'-','color',mycolor(1,:)); hold on
    end
    
    mytitle = ['$q_3$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = ''; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot Q4
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,q(:,4),'.','markersize',14,'color',mycolor(1,:)); hold on
    
    plot(Wie.tArray,Wie.qArray(4,:),'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,q4l,'-','color',mycolor(1,:)); hold on
    end
    
    mytitle = ['$q_4$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = ''; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot w1
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,w(:,1),'.','markersize',14,'color',mycolor(6,:)); hold on
    
    plot(Wie.tArray,Wie.wArray(1,:),'color',mycolor(2,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,w1l,'-','color',mycolor(6,:)); hold on
    end
    
    mytitle = ['$\omega_1$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'rad/sec'; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot w2
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,w(:,2),'.','markersize',14,'color',mycolor(6,:)); hold on
    
    plot(Wie.tArray,Wie.wArray(2,:),'color',mycolor(2,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,w2l,'-','color',mycolor(6,:)); hold on
    end
    
    mytitle = ['$\omega_2$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'rad/sec'; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot w3
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longt,w(:,3),'.','markersize',14,'color',mycolor(6,:)); hold on
    
    plot(Wie.tArray,Wie.wArray(3,:),'color',mycolor(2,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,w3l,'-','color',mycolor(6,:)); hold on
    end
    
    mytitle = ['$\omega_3$']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'rad/sec'; % y label with latex
    mylegend = {'Gauss PS Method','SR Robust'}; % legend with latex
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
    % Plot gimbal angles
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
    
    % Convert to degrees
    delta = delta*180/pi;
    d1l = d1l*180/pi;
    d2l = d2l*180/pi;
    d3l = d3l*180/pi;
    
    if strcmp(method,'Pseudospectral')
        plot(T,d1l,'-','color',mycolor(1,:)); hold on
        plot(T,d2l,':','color',mycolor(4,:),'linewidth',2); hold on
        plot(T,d3l,'--','color',mycolor(2,:)); hold on
    end
        
    plot(longt,delta(:,1),'.','markersize',14,'color',mycolor(1,:)); hold on
    plot(longt,delta(:,2),'.','markersize',14,'color',mycolor(4,:)); hold on
    plot(longt,delta(:,3),'.','markersize',14,'color',mycolor(2,:)); hold on
    
    mytitle = ['Gimbal Angles']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'degrees'; % y label with latex
    mylegend = {'$\delta_1$','$\delta_2$','$\delta_3$'}; % legend with latex
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
    % Plot Jacobian determinant
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
    
    plot(longt,Det,'.','markersize',14,'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,Detl,'-','color',mycolor(7,:)); hold on
    end
    
    mytitle = ['Determinant of Jacobian']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = '$\det(A)$'; % y label with latex
%     myylabel = '';
    xlabel(myxlabel) % create x label
    ylabel(myylabel) % create y label
    ylim([-1.5e9,1.5e9])
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
        
    
end