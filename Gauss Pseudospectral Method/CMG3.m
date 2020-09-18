%--------------------------------------------------------------------------
% CMG3.m
% Solve the 3/4-CMG pyramid configuration using
% the differential Gauss pseudospectral method
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

% Run Wie simulation, for comparison
WieSimulation

% Problem parameters
p.ns = 10; p.nu = 3; % number of states and controls
p.t0 = 0;            % time horizon (final time to be determined)
p.j = diag([21400,20100,5000]); % spacecraft tensor of inertia, in SI units
p.jInv = inv(p.j);
p.beta = 53.13; % skew angle, in degrees
p.h = 1000; % each CMG momentum, in SI units (constant)
slewMax = 10; % max slew rate, in deg/s
p.slewMax = slewMax*pi/180; % convert to rad/s
p.deltaMax = 2; % max gimbal rate, in rad/s
p.q0 = [0; 0; 0; 1]; % quaternion initial condition
p.qf = [0.4; 0; 0; sqrt(0.84)]; % quaternion final condition
p.w0 = [0; 0; 0]; % angular velocity initial condition
p.wf = [0; 0; 0]; % angular velocity final condition
delta0 = [60; 180; -60]; % initial gimbal angles, in degrees
% delta0 = [90; 0; -90]; % singularity
p.delta0 = delta0*pi/180; % convert to radians
p.tfGuess = 10; % guess of final time, in seconds
p.tfUpper = 12; % upper bound on maneuver time

% Direct transcription parameters
% pre-computed values for tau,w,D available for N = 10,20,50,100
% p.nt = 10; % number of node points, including -1
p.nt = 20; % number of node points, including -1
% p.nt = 50; % number of node points, including -1
p.tau = double(Gauss_nodes(p.nt-1)); % scaled time horizon
p.D =  double(Gauss_Dmatrix(p.tau)); % differential approximation matrix
p.W = double(Gauss_weights(p.tau)); % for Gaussian quadrature

% Discretized variable indices in x = [q1;...w1;...d1;...;u1...];
n = p.nt;
p.q1i = 1:n; p.q2i = n+1:2*n;
p.q3i = 2*n+1:3*n; p.q4i = 3*n+1:4*n; m = 4*n;
p.w1i = m+1:m+n; p.w2i = m+n+1:m+2*n;
p.w3i = m+2*n+1:m+3*n; m = m+3*n;
p.d1i = m+1:m+n; p.d2i = m+n+1:m+2*n;
p.d3i = m+2*n+1:m+3*n; m = m+3*n;
p.u1i = m+1:m+n-1; p.u2i = m+n:m+2*n-2;
p.u3i = m+2*n-1:m+3*n-3;
p.tfi = m+3*n-2;
len = m+3*n-2; % length of x

% Initial guess
x0 = initialGuess(Wie,p,len);

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6, ...
    'Algorithm','sqp','MaxIterations',1000);

% Solve the problem
% x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p,2),options);
% x0 = x;
% x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p,3),options);
% x0 = x;
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p,1),options);

% temporary
% x = x0;

% Extract the optimal solution
q1 = x(p.q1i); q2 = x(p.q2i); q3 = x(p.q3i); q4 = x(p.q4i);
w1 = x(p.w1i); w2 = x(p.w2i); w3 = x(p.w3i);
d1 = x(p.d1i); d2 = x(p.d2i); d3 = x(p.d3i);
u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
q = [q1, q2, q3, q4];
w = [w1, w2, w3];
delta = [d1, d2, d3];
u = [u1, u2, u3];
p.tf = x(p.tfi); % extract final time
p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon

% Calculate determinant of Jacobian over time
Det = zeros(length(d1),1);
for it = 1:length(d1)
    Det(it) = det(Jacobian3(p.h, p.beta, [d1(it);d2(it);d3(it)]*180/pi));
end

% Obtain the costate
% Lam = zeros(p.nt + 1, 1);
% Lam(1) = lambda.eqnonlin(1);
% Lam(end) = -lambda.eqnonlin(2);
% Lam(2:end-1) = lambda.eqnonlin(3:end)./p.w + Lam(end);
% Lam = -Lam;

% Plots
% Plots_Gauss_CMG3(q,w,delta,u,Det,p,'Pseudospectral')
Plots2_Gauss_CMG3(q,w,delta,u,Det,p,Wie,'Pseudospectral') % for paper

% Initial guess
function x0 = initialGuess(Wie,p,len)
    x0 = zeros(len,1);
    indices = round(0.5*(p.tau+1)*p.tfGuess/Wie.dt);
    indices = [1; indices];
    x0(p.q1i) = Wie.qArray(1,indices)';
    x0(p.q2i) = Wie.qArray(2,indices)';
    x0(p.q3i) = Wie.qArray(3,indices)';
    x0(p.q4i) = Wie.qArray(4,indices)';
    x0(p.w1i) = Wie.wArray(1,indices)';
    x0(p.w2i) = Wie.wArray(2,indices)';
    x0(p.w3i) = Wie.wArray(3,indices)';
    x0(p.tfi) = p.tfGuess;
    x0(p.d1i) = p.delta0(1);
    x0(p.d2i) = p.delta0(2);
    x0(p.d3i) = p.delta0(3);
    
    % temporary
%     x0(p.w1i) = x0(p.w1i)*180/pi;
%     x0(p.w2i) = x0(p.w2i)*180/pi;
%     x0(p.w3i) = x0(p.w3i)*180/pi;
end

% Objective function
function J = objective(x,p)
    % Extract
    tf = x(p.tfi); % extract final time
    d1 = x(p.d1i); d2 = x(p.d2i); d3 = x(p.d3i); % extract gimbal angles
    
    % Calculate determinant of Jacobian over time
    Det = zeros(length(d1),1);
    for it = 1:length(d1)
        Det(it) = det(Jacobian3(p.h, p.beta, [d1(it);d2(it);d3(it)]*180/pi));
    end
    Det = Det(2:end);
    
    % Weights
    Weight1 = 0; % for final time
%     Weight2 = 0; % for maximal integrated det(A)
    Weight2 = 1e-19;
    
    L = -Weight2*abs(Det); % integrand
    J = Weight1*tf + ((tf-p.t0)/2)*dot(p.W,L); % cost function
end

% Constraint function
function [c,ceq] = constraints(x,p,choice)
    % Extract
    q1 = x(p.q1i); q2 = x(p.q2i); q3 = x(p.q3i); q4 = x(p.q4i);
    w1 = x(p.w1i); w2 = x(p.w2i); w3 = x(p.w3i);
    d1 = x(p.d1i); d2 = x(p.d2i); d3 = x(p.d3i);
    u1 = x(p.u1i); u2 = x(p.u2i); u3 = x(p.u3i);
    tf = x(p.tfi); % extract final time
    
    qArray = [q1'; q2'; q3'; q4'];
    wArray = [w1'; w2'; w3'];
    deltaArray = [d1'; d2'; d3'];
    uArray = [u1'; u2'; u3'];
    
    qvDotArray = zeros(3, p.nt-1);
    q4DotArray = zeros(1, p.nt-1);
    wDotArray = zeros(3, p.nt-1);
    
    % Equations of motion
    for it = 2:p.nt
        qvDotArray(:,it-1) = -0.5*cross(wArray(:,it), qArray(1:3,it)) ...
            + 0.5*qArray(4,it)*wArray(:,it);
        q4DotArray(it-1) = -0.5*dot(wArray(:,it), qArray(1:3,it));
        
        hv = angVector3(p.h, p.beta, deltaArray(:,it)*180/pi);
        A = Jacobian3(p.h, p.beta, deltaArray(:,it)*180/pi);
        
        wDotArray(:,it-1) = -p.jInv*(cross(wArray(:,it), p.j*wArray(:,it)+hv) ...
            + A*uArray(:,it-1));
    end
    
    % Create matrices: (p.nt x p.ns) for Y, one fewer row for F.
    % For Gauss pseudospectral, k=0 is needed here:
    Y = [q1,q2,q3,q4,w1,w2,w3,d1,d2,d3];
    % For Gauss PS, no k=0 here:
    F = ((tf-p.t0)/2)*[qvDotArray', q4DotArray', wDotArray', uArray'];
    
    % Initial state conditions
    ceq1 = qArray(:,1) - p.q0;
    ceq2 = wArray(:,1) - p.w0;
    ceq3 = deltaArray(:,1) - p.delta0;

    % Final state conditions
    % necessary for *Gauss* pseudospectral method
    ceq4 = qArray(:,1) + dot([p.W,p.W,p.W,p.W],F(:,1:4))' - p.qf;
    ceq5 = wArray(:,1) + dot([p.W,p.W,p.W],F(:,5:7))' - p.wf;
    
    ceq6 = p.D*Y - F; % differential equation approximation
    
    ceq7 = vecnorm(qArray)' - 1; % quaternion condition
    
    % Path constraints
    c1 = abs(u1) - p.deltaMax;
    c2 = abs(u2) - p.deltaMax;
    c3 = abs(u3) - p.deltaMax;
    c4 = vecnorm(wArray)' - p.slewMax;
    c5 = -tf;
    c6 = tf - p.tfUpper;
    
    % Combine constraints
    if choice == 0
        c = [c1;c2;c3;c4;c5;c6];
        ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6(:);ceq7];
    elseif choice == 1
        c = [c1;c2;c3;c4;c5;c6];
        ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6(:)];
    elseif choice == 2
        c = [c1;c2;c3;c4;c5;c6];
        ceq = [ceq1;ceq2;ceq3;ceq4;ceq6(:)];
    elseif choice == 3
        c = [c1;c2;c3;c4;c5;c6];
        ceq = [ceq1;ceq2;ceq3;ceq4;ceq6(:);ceq7];
    end
end