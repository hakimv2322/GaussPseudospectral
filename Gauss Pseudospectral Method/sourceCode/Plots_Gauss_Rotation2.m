%--------------------------------------------------------------------------
% Plots_Gauss_Rotation2.m
% Plotting function, Gauss pseudospectral method
% For rotation problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Underlying code: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
function Plots_Gauss_Rotation2(qA,wA,uA,qB,wB,uB,Lam,p,method)
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
    
    LamqA = Lam.LamqA; LamwA = Lam.LamwA;
    LamqB = Lam.LamqB; LamwB = Lam.LamwB;
    muA = Lam.muA; muB = Lam.muB;
    
    if strcmp(method,'Pseudospectral')
        q1lA = LagrangeInter(longtA',qA(:,1)',TA);
        q2lA = LagrangeInter(longtA',qA(:,2)',TA);
        q3lA = LagrangeInter(longtA',qA(:,3)',TA);
        q4lA = LagrangeInter(longtA',qA(:,4)',TA);
        w1lA = LagrangeInter(longtA',wA(:,1)',TA);
        w2lA = LagrangeInter(longtA',wA(:,2)',TA);
        w3lA = LagrangeInter(longtA',wA(:,3)',TA);
        u1lA = LagrangeInter(p.tA',uA(:,1)',TA);
        u2lA = LagrangeInter(p.tA',uA(:,2)',TA);
        u3lA = LagrangeInter(p.tA',uA(:,3)',TA);
        Lamq1lA = LagrangeInter(longertA',LamqA(1,:),TA);
        Lamq2lA = LagrangeInter(longertA',LamqA(2,:),TA);
        Lamq3lA = LagrangeInter(longertA',LamqA(3,:),TA);
        Lamq4lA = LagrangeInter(longertA',LamqA(4,:),TA);
        Lamw1lA = LagrangeInter(longertA',LamwA(1,:),TA);
        Lamw2lA = LagrangeInter(longertA',LamwA(2,:),TA);
        Lamw3lA = LagrangeInter(longertA',LamwA(3,:),TA);
        mu1lA = LagrangeInter(p.tA',muA(1,:),TA);
        mu2lA = LagrangeInter(p.tA',muA(2,:),TA);
        mu3lA = LagrangeInter(p.tA',muA(3,:),TA);
        q1lB = LagrangeInter(longtB',qB(:,1)',TB);
        q2lB = LagrangeInter(longtB',qB(:,2)',TB);
        q3lB = LagrangeInter(longtB',qB(:,3)',TB);
        q4lB = LagrangeInter(longtB',qB(:,4)',TB);
        w1lB = LagrangeInter(longtB',wB(:,1)',TB);
        w2lB = LagrangeInter(longtB',wB(:,2)',TB);
        w3lB = LagrangeInter(longtB',wB(:,3)',TB);
        u1lB = LagrangeInter(p.tB',uB(:,1)',TB);
        u2lB = LagrangeInter(p.tB',uB(:,2)',TB);
        u3lB = LagrangeInter(p.tB',uB(:,3)',TB);
        Lamq1lB = LagrangeInter(longertB',LamqB(1,:),TB);
        Lamq2lB = LagrangeInter(longertB',LamqB(2,:),TB);
        Lamq3lB = LagrangeInter(longertB',LamqB(3,:),TB);
        Lamq4lB = LagrangeInter(longertB',LamqB(4,:),TB);
        Lamw1lB = LagrangeInter(longertB',LamwB(1,:),TB);
        Lamw2lB = LagrangeInter(longertB',LamwB(2,:),TB);
        Lamw3lB = LagrangeInter(longertB',LamwB(3,:),TB);
        mu1lB = LagrangeInter(p.tB',muB(1,:),TB);
        mu2lB = LagrangeInter(p.tB',muB(2,:),TB);
        mu3lB = LagrangeInter(p.tB',muB(3,:),TB);
        q1l = [q1lA, q1lB(2:end)];
        q2l = [q2lA, q2lB(2:end)];
        q3l = [q3lA, q3lB(2:end)];
        q4l = [q4lA, q4lB(2:end)];
        w1l = [w1lA, w1lB(2:end)];
        w2l = [w2lA, w2lB(2:end)];
        w3l = [w3lA, w3lB(2:end)];
        u1l = [u1lA, u1lB(2:end)];
        u2l = [u2lA, u2lB(2:end)];
        u3l = [u3lA, u3lB(2:end)];
        Lamq1l = [Lamq1lA, Lamq1lB(2:end)];
        Lamq2l = [Lamq2lA, Lamq2lB(2:end)];
        Lamq3l = [Lamq3lA, Lamq3lB(2:end)];
        Lamq4l = [Lamq4lA, Lamq4lB(2:end)];
        Lamw1l = [Lamw1lA, Lamw1lB(2:end)];
        Lamw2l = [Lamw2lA, Lamw2lB(2:end)];
        Lamw3l = [Lamw3lA, Lamw3lB(2:end)];
        mu1l = [mu1lA, mu1lB(2:end)];
        mu2l = [mu2lA, mu2lB(2:end)];
        mu3l = [mu3lA, mu3lB(2:end)];
    end
    
    q = [qA; qB];
    w = [wA; wB];
    u = [uA; uB];
    Lamq = [LamqA'; LamqB(:,2:end)'];
    Lamw = [LamwA'; LamwB(:,2:end)'];
    
    
    %----------------------------------------------------------------------
    % plot quaternion magnitude
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
    
    qMag = sqrt(q1l.^2 + q2l.^2 + q3l.^2 + q4l.^2);
    plot(T,qMag,'-','color',mycolor(7,:),'linewidth',1); hold on
    
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
    % plot quaternion costates
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longert,Lamq(:,1),'.','markersize',14,'color',mycolor(1,:)); hold on
    plot(longert,Lamq(:,2),'.','markersize',14,'color',mycolor(3,:)); hold on
    plot(longert,Lamq(:,3),'.','markersize',14,'color',mycolor(5,:)); hold on
    plot(longert,Lamq(:,4),'.','markersize',14,'color',mycolor(7,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,Lamq1l,'-','color',mycolor(1,:)); hold on
        plot(T,Lamq2l,'-','color',mycolor(3,:)); hold on
        plot(T,Lamq3l,'-','color',mycolor(5,:)); hold on
        plot(T,Lamq4l,'-','color',mycolor(7,:)); hold on
    end
    
    mytitle = ['Quaternion Costates']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'Quaternion Costates'; % y label with latex
    mylegend = {'$\Lambda$q1','$\Lambda$q2','$\Lambda$q3','$\Lambda$q4'}; % legend with latex
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
    % plot angular velocity costates
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color
        
    plot(longert,Lamw(:,1),'.','markersize',14,'color',mycolor(2,:)); hold on
    plot(longert,Lamw(:,2),'.','markersize',14,'color',mycolor(4,:)); hold on
    plot(longert,Lamw(:,3),'.','markersize',14,'color',mycolor(6,:)); hold on
    
    if strcmp(method,'Pseudospectral')
        plot(T,Lamw1l,'-','color',mycolor(2,:)); hold on
        plot(T,Lamw2l,'-','color',mycolor(4,:)); hold on
        plot(T,Lamw3l,'-','color',mycolor(6,:)); hold on
    end
    
    mytitle = ['Angular Velocity Costates']; % title with latex
    myxlabel = '$t$ (sec.)'; % x label with latex
    myylabel = 'Angular Velocity Costates'; % y label with latex
    mylegend = {'$\Lambda\omega_x$','$\Lambda\omega_y$','$\Lambda\omega_z$'}; % legend with latex
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
    
    %----------------------------------------------------------------------
    % Calculate and plot Hamiltonian
    f = calculateF(q1l,q2l,q3l,q4l,w1l,w2l,w3l,u1l,u2l,u3l,p);
    Laml = [Lamq1l;Lamq2l;Lamq3l;Lamq4l;Lamw1l;Lamw2l;Lamw3l];
    mul = [mu1l;mu2l;mu3l];
    Cl = [abs(u1l) - p.u1max; abs(u2l) - p.u2max; abs(u3l) - p.u3max;];
    Ham = 1 + dot(Laml,f) - dot(mul,Cl);
    
    hf = figure; % create a new figure and save handle
    hf.Color = wcolor; % change the figure background color

    if strcmp(method,'Pseudospectral')
        plot(T,Ham,'-','color',mycolor(1,:),'linewidth',1); hold on
    end
    
	mytitle = ['Gauss ',method,' - Hamiltonian']; % title with latex
    myxlabel = '$t$'; % x label with latex
    myylabel = 'Hamiltonian'; % y label with latex
%     ylim([-1e8, 1e8])
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

% Calculate F for the Hamiltonian
function F = calculateF(q1,q2,q3,q4,w1,w2,w3,u1,u2,u3,p)
    % assumes all inputs are row vectors, so convert into columns
    q1 = q1'; q2 = q2'; q3 = q3'; q4 = q4';
    w1 = w1'; w2 = w2'; w3 = w3';
    u1 = u1'; u2 = u2'; u3 = u3';
    
    w = [w1,w2,w3];
    q = [q1,q2,q3,q4];
    u = [u1,u2,u3];
    zer = zeros(length(w1),1);
    
    % Omega matrix
    Om1 = [zer w3 -w2 w1];
    Om2 = [-w3 zer w1 w2];
    Om3 = [w2 -w1 zer w3];
    Om4 = [-w1 -w2 -w3 zer];
    
    % Dynamic equations
    dq1 = diag(0.5*Om1*q');
    dq2 = diag(0.5*Om2*q');
    dq3 = diag(0.5*Om3*q');
    dq4 = diag(0.5*Om4*q');
    
    dw = p.jInv*(u' - cross(w',p.j*w'));
    dw = dw';
    
    F = [dq1,dq2,dq3,dq4,dw]';
end
    
  