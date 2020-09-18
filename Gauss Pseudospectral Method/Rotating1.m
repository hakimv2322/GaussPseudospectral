%--------------------------------------------------------------------------
% Rotating1.m
% Solve a rotations problem using the differential Gauss
% pseudospectral method
% (namely Gauss nodes and Gaussian quadrature)
% Method developed by David Benson, MIT PhD thesis
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
% Underlying source code by:
% Daniel R. Herber, Graduate Student, University
% of Illinois at Urbana-Champaign, website:
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------

close all
clear all

% problem parameters
p.ns = 7; p.nu = 3; % number of states and controls
p.t0 = 0;           % time horizon (final time to be determined)
p.j = diag([150,100,40]); % Tensor of inertia for rigid body
p.jInv = inv(p.j);
p.u1max = 1; % Maximum control torque
p.u2max = 1; % Maximum control torque
p.u3max = 1; % Maximum control torque
q0 = [0 0 0 1]; % Quaternion initial condition
qf = 0.5*sqrt(2)*[0 0 1 1]; % Quaternion final condition
% qf = [1/sqrt(6) 1/sqrt(6) 1/sqrt(6) 1/sqrt(2)]; % Quaternion final condition
w0 = [0 0 0]; % Angular velocity initial condition
wf = [0 0 0]; % Angular velocity final condition
% States 1-4 are the quaternions, 5-7 are the angular velocity:
p.y10 = q0(1); p.y20 = q0(2); p.y30 = q0(3); p.y40 = q0(4);
p.y1f = qf(1); p.y2f = qf(2); p.y3f = qf(3); p.y4f = qf(4);
p.y50 = w0(1); p.y60 = w0(2); p.y70 = w0(3);
p.y5f = wf(1); p.y6f = wf(2); p.y7f = wf(3);

% direct transcription parameters
% pre-computed values for tau,w,D available for N = 10,20,50,100
p.nt = 10; % number of node points, including -1
% p.nt = 20; % number of node points; this takes a few minutes
% p.nt = 50; % number of node points; this takes a few minutes
p.tau = double(Gauss_nodes(p.nt-1)); % scaled time horizon
p.D =  double(Gauss_Dmatrix(p.tau)); % differential approximation matrix
p.w = double(Gauss_weights(p.tau)); % for gaussian quadrature
% discretized variable indices in x = [y1,y2,...,y7,u1,u2,u3];
n = p.nt;
p.y1i = 1:n; p.y2i = n+1:2*n; p.y3i = 2*n+1:3*n; p.y4i = 3*n+1:4*n;
p.y5i = 4*n+1:5*n; p.y6i = 5*n+1:6*n; p.y7i = 6*n+1:7*n;
p.u1i = 7*n+1:8*n-1; p.u2i = 8*n:9*n-2; p.u3i = 9*n-1:10*n-3;
p.tfi = 10*n-2; % Be sure to add 1 in the zeros below:
x0 = zeros(p.nt*(p.ns+p.nu) - p.nu + 1,1); % initial guess (all zeros)
x0(p.tfi) = 1;
x0(p.y1i) = p.y10; x0(p.y2i) = p.y20; x0(p.y3i) = p.y30; x0(p.y4i) = p.y40;
x0(p.y5i) = p.y50; x0(p.y6i) = p.y60; x0(p.y7i) = p.y70;

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6,'Algorithm','sqp');

% Find a "feasible point" for the initial guess
x = fmincon(@(x) 1,x0,[],[],[],[],[],[],@(x) constraintsWeird(x,p),options);
x0(p.y1i) = x(p.y1i); x0(p.y2i) = x(p.y2i); x0(p.y3i) = x(p.y3i); x0(p.y4i) = x(p.y4i);
x0(p.y5i) = x(p.y5i); x0(p.y6i) = x(p.y6i); x0(p.y7i) = x(p.y7i);
x0(p.tfi) = x(p.tfi);

% Find a "feasible point" for the initial guess (iteration 1.5)
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraintsWeird(x,p),options);
x0(p.y1i) = x(p.y1i); x0(p.y2i) = x(p.y2i); x0(p.y3i) = x(p.y3i); x0(p.y4i) = x(p.y4i);
x0(p.y5i) = x(p.y5i); x0(p.y6i) = x(p.y6i); x0(p.y7i) = x(p.y7i);
x0(p.tfi) = x(p.tfi);

% Find a "feasible point" for the initial guess (2nd iteration)
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraintsWeird2(x,p),options);
x0(p.y1i) = x(p.y1i); x0(p.y2i) = x(p.y2i); x0(p.y3i) = x(p.y3i); x0(p.y4i) = x(p.y4i);
x0(p.y5i) = x(p.y5i); x0(p.y6i) = x(p.y6i); x0(p.y7i) = x(p.y7i);
x0(p.tfi) = x(p.tfi);

% Find a "feasible point" for the initial guess (3rd iteration)
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraintsWeird3(x,p),options);
x0(p.y1i) = x(p.y1i); x0(p.y2i) = x(p.y2i); x0(p.y3i) = x(p.y3i); x0(p.y4i) = x(p.y4i);
x0(p.y5i) = x(p.y5i); x0(p.y6i) = x(p.y6i); x0(p.y7i) = x(p.y7i);
x0(p.tfi) = x(p.tfi);

% Find a "feasible point" for the initial guess (4th iteration)
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraintsWeird4(x,p),options);
x0(p.y1i) = x(p.y1i); x0(p.y2i) = x(p.y2i); x0(p.y3i) = x(p.y3i); x0(p.y4i) = x(p.y4i);
zer = zeros(length(p.y5i),1);
x0(p.y5i) = zer; x0(p.y6i) = zer; x0(p.y7i) = x(p.y7i);
x0(p.tfi) = x(p.tfi);

% solve the problem; obtain the optimal solution
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);
% Extract:
y1 = x(p.y1i); y2 = x(p.y2i); y3 = x(p.y3i); y4 = x(p.y4i);
y5 = x(p.y5i); y6 = x(p.y6i); y7 = x(p.y7i);
u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
q = [y1, y2, y3, y4];
w = [y5, y6, y7];
u = [u1, u2, u3];
p.tf = x(p.tfi); % extract final time
p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon
% plots
Plots_Gauss_Rotation1(q,w,u,p,'Pseudospectral')

% objective function
function J = objective(x,p)
    % Extract:
    tf = x(p.tfi);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    uSq = u1.^2 + u2.^2 + u3.^2;
    
    % Calculate objective:
    J = tf; % for minimum time
%     J = ((tf-p.t0)/2)*dot(p.w,uSq); % for minimum control effort
end

% constraint function
function [c,ceq] = constraints(x,p)
    % Extract:
    q1 = x(p.y1i); q2 = x(p.y2i); q3 = x(p.y3i); q4 = x(p.y4i);
    w1 = x(p.y5i); w2 = x(p.y6i); w3 = x(p.y7i);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    tf = x(p.tfi); % extract final time
    
    w = [w1(2:end),w2(2:end),w3(2:end)];
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
    
    % create matrices: (p.nt x p.ns) for Y, one fewer row for F.
    % For Gauss pseudospectral, k=0 is needed here:
    Y = [q1,q2,q3,q4,w1,w2,w3];
    % For Gauss PS, no k=0 here:
    F0 = [dq1(2:end),dq2(2:end),dq3(2:end),dq4(2:end),dw];
    F = ((tf-p.t0)/2)*F0;
    
    % initial state conditions
    ceq1 = q1(1) - p.y10;
    ceq2 = q2(1) - p.y20;
    ceq3 = q3(1) - p.y30;
    ceq4 = q4(1) - p.y40;
    ceq5 = w1(1) - p.y50;
    ceq6 = w2(1) - p.y60;
    ceq7 = w3(1) - p.y70;
    
    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq8 = q1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq9 = q2(1) + dot(p.w,F(:,2)) - p.y2f;
    ceq10 = q3(1) + dot(p.w,F(:,3)) - p.y3f;
    ceq11 = q4(1) + dot(p.w,F(:,4)) - p.y4f;
    ceq12 = w1(1) + dot(p.w,F(:,5)) - p.y5f;
    ceq13 = w2(1) + dot(p.w,F(:,6)) - p.y6f;
    ceq14 = w3(1) + dot(p.w,F(:,7)) - p.y7f;

    ceq15 = p.D*Y - F; % differential equation approximation
    ceq16 = sqrt(q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1; % quaternion condition
    
    % path constraints
    c1 = abs(u1) - p.u1max;
    c2 = abs(u2) - p.u2max;
    c3 = abs(u3) - p.u3max;
    c4 = -tf;
    
    % combine constraints
    c = [c1;c2;c3;c4;];
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7;ceq8; ...
        ceq9;ceq10;ceq11;ceq12;ceq13;ceq14;ceq15(:);];
end


% constraint function (weird), for finding a feasible initial guess
function [c,ceq] = constraintsWeird(x,p)
    % Extract:
    q1 = x(p.y1i); q2 = x(p.y2i); q3 = x(p.y3i); q4 = x(p.y4i);
    w1 = x(p.y5i); w2 = x(p.y6i); w3 = x(p.y7i);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    tf = x(p.tfi); % extract final time
    
    w = [w1(2:end),w2(2:end),w3(2:end)];
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
    
    % create matrices (p.nt x p.ns)
    % For Gauss pseudospectral, k=0 is needed here:
    Y = [q1,q2,q3,q4,w1,w2,w3];
    % For Gauss PS, no k=0 here:
    F0 = [dq1(2:end),dq2(2:end),dq3(2:end),dq4(2:end),dw'];
    F = ((tf-p.t0)/2)*F0;
    
    % initial state conditions
    ceq1 = q1(1) - p.y10;
    ceq2 = q2(1) - p.y20;
    ceq3 = q3(1) - p.y30;
    ceq4 = q4(1) - p.y40;
    ceq5 = w1(1) - p.y50;
    ceq6 = w2(1) - p.y60;
    ceq7 = w3(1) - p.y70;
    
    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq8 = q1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq9 = q2(1) + dot(p.w,F(:,2)) - p.y2f;
    ceq10 = q3(1) + dot(p.w,F(:,3)) - p.y3f;
    ceq11 = q4(1) + dot(p.w,F(:,4)) - p.y4f;
    ceq12 = w1(1) + dot(p.w,F(:,5)) - p.y5f;
    ceq13 = w2(1) + dot(p.w,F(:,6)) - p.y6f;
    ceq14 = w3(1) + dot(p.w,F(:,7)) - p.y7f;

    ceq15 = p.D*Y - F; % differential equation approximation
    ceq16 = sqrt(q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1; % quaternion condition
    
    % path constraints
    c1 = abs(u1) - p.u1max;
    c2 = abs(u2) - p.u2max;
    c3 = abs(u3) - p.u3max;
    c4 = -tf;
    
    % combine constraints
    c = [c1;c2;c3;c4;];
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq8; ...
        ceq9;ceq10;ceq11;ceq15(:);ceq16;];
end

% constraint function (weird2), for finding a feasible initial guess
function [c,ceq] = constraintsWeird2(x,p)
    % Extract:
    q1 = x(p.y1i); q2 = x(p.y2i); q3 = x(p.y3i); q4 = x(p.y4i);
    w1 = x(p.y5i); w2 = x(p.y6i); w3 = x(p.y7i);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    tf = x(p.tfi); % extract final time
    
    w = [w1(2:end),w2(2:end),w3(2:end)];
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
    
    % create matrices: (p.nt x p.ns) for Y, one fewer row for F.
    % For Gauss pseudospectral, k=0 is needed here:
    Y = [q1,q2,q3,q4,w1,w2,w3];
    % For Gauss PS, no k=0 here:
    F0 = [dq1(2:end),dq2(2:end),dq3(2:end),dq4(2:end),dw'];
    F = ((tf-p.t0)/2)*F0;
    
    % initial state conditions
    ceq1 = q1(1) - p.y10;
    ceq2 = q2(1) - p.y20;
    ceq3 = q3(1) - p.y30;
    ceq4 = q4(1) - p.y40;
    ceq5 = w1(1) - p.y50;
    ceq6 = w2(1) - p.y60;
    ceq7 = w3(1) - p.y70;
    
    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq8 = q1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq9 = q2(1) + dot(p.w,F(:,2)) - p.y2f;
    ceq10 = q3(1) + dot(p.w,F(:,3)) - p.y3f;
    ceq11 = q4(1) + dot(p.w,F(:,4)) - p.y4f;
    ceq12 = w1(1) + dot(p.w,F(:,5)) - p.y5f;
    ceq13 = w2(1) + dot(p.w,F(:,6)) - p.y6f;
    ceq14 = w3(1) + dot(p.w,F(:,7)) - p.y7f;

    ceq15 = p.D*Y - F; % differential equation approximation
    ceq16 = sqrt(q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1; % quaternion condition
    
    % path constraints
    c1 = abs(u1) - p.u1max;
    c2 = abs(u2) - p.u2max;
    c3 = abs(u3) - p.u3max;
    c4 = -tf;
    
    % combine constraints
    c = [c1;c2;c3;c4;];
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq8; ...
        ceq9;ceq10;ceq11;ceq12;ceq13;ceq15(:);ceq16;];
end

% constraint function (weird3), for finding a feasible initial guess
function [c,ceq] = constraintsWeird3(x,p)
    % Extract:
    q1 = x(p.y1i); q2 = x(p.y2i); q3 = x(p.y3i); q4 = x(p.y4i);
    w1 = x(p.y5i); w2 = x(p.y6i); w3 = x(p.y7i);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    tf = x(p.tfi); % extract final time
    
    w = [w1(2:end),w2(2:end),w3(2:end)];
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
    
    % create matrices: (p.nt x p.ns) for Y, one fewer row for F.
    % For Gauss pseudospectral, k=0 is needed here:
    Y = [q1,q2,q3,q4,w1,w2,w3];
    % For Gauss PS, no k=0 here:
    F0 = [dq1(2:end),dq2(2:end),dq3(2:end),dq4(2:end),dw'];
    F = ((tf-p.t0)/2)*F0;
    
    % initial state conditions
    ceq1 = q1(1) - p.y10;
    ceq2 = q2(1) - p.y20;
    ceq3 = q3(1) - p.y30;
    ceq4 = q4(1) - p.y40;
    ceq5 = w1(1) - p.y50;
    ceq6 = w2(1) - p.y60;
    ceq7 = w3(1) - p.y70;
    
    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq8 = q1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq9 = q2(1) + dot(p.w,F(:,2)) - p.y2f;
    ceq10 = q3(1) + dot(p.w,F(:,3)) - p.y3f;
    ceq11 = q4(1) + dot(p.w,F(:,4)) - p.y4f;
    ceq12 = w1(1) + dot(p.w,F(:,5)) - p.y5f;
    ceq13 = w2(1) + dot(p.w,F(:,6)) - p.y6f;
    ceq14 = w3(1) + dot(p.w,F(:,7)) - p.y7f;

    ceq15 = p.D*Y - F; % differential equation approximation
    ceq16 = sqrt(q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1; % quaternion condition
    
    % path constraints
    c1 = abs(u1) - p.u1max;
    c2 = abs(u2) - p.u2max;
    c3 = abs(u3) - p.u3max;
    c4 = -tf;
    
    % combine constraints
    c = [c1;c2;c3;c4;];
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7;ceq8; ...
        ceq9;ceq10;ceq11;ceq12;ceq13;ceq15(:);];
end

% constraint function (weird4), for finding a feasible initial guess
function [c,ceq] = constraintsWeird4(x,p)
    % Extract:
    q1 = x(p.y1i); q2 = x(p.y2i); q3 = x(p.y3i); q4 = x(p.y4i);
    w1 = x(p.y5i); w2 = x(p.y6i); w3 = x(p.y7i);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    tf = x(p.tfi); % extract final time
    
    w = [w1(2:end),w2(2:end),w3(2:end)];
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
    
    % create matrices: (p.nt x p.ns) for Y, one fewer row for F.
    % For Gauss pseudospectral, k=0 is needed here:
    Y = [q1,q2,q3,q4,w1,w2,w3];
    % For Gauss PS, no k=0 here:
    F0 = [dq1(2:end),dq2(2:end),dq3(2:end),dq4(2:end),dw'];
    F = ((tf-p.t0)/2)*F0;
    
    % initial state conditions
    ceq1 = q1(1) - p.y10;
    ceq2 = q2(1) - p.y20;
    ceq3 = q3(1) - p.y30;
    ceq4 = q4(1) - p.y40;
    ceq5 = w1(1) - p.y50;
    ceq6 = w2(1) - p.y60;
    ceq7 = w3(1) - p.y70;
    
    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq8 = q1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq9 = q2(1) + dot(p.w,F(:,2)) - p.y2f;
    ceq10 = q3(1) + dot(p.w,F(:,3)) - p.y3f;
    ceq11 = q4(1) + dot(p.w,F(:,4)) - p.y4f;
    ceq12 = w1(1) + dot(p.w,F(:,5)) - p.y5f;
    ceq13 = w2(1) + dot(p.w,F(:,6)) - p.y6f;
    ceq14 = w3(1) + dot(p.w,F(:,7)) - p.y7f;

    ceq15 = p.D*Y - F; % differential equation approximation
    ceq16 = sqrt(q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1; % quaternion condition
    
    % path constraints
    c1 = abs(u1) - p.u1max;
    c2 = abs(u2) - p.u2max;
    c3 = abs(u3) - p.u3max;
    c4 = -tf;
    
    % combine constraints
    c = [c1;c2;c3;c4;];
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7;ceq8; ...
        ceq9;ceq10;ceq11;ceq12;ceq13;ceq14;ceq15(:);];
end

