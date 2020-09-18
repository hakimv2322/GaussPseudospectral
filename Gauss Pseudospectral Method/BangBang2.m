%--------------------------------------------------------------------------
% BangBang2.m
% (Using two-phase approach)
% Solve the Bang-Bang problem using the differential Gauss
% pseudospectral method, p. 141 of Benson's PhD thesis
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
p.ns = 2; p.nu = 1; % number of states and controls
p.t0 = 0;           % time horizon (final time to be determined)
p.y10 = 1; p.y1f = 0; p.y20 = 3; p.y2f = 0; % boundary conditions

% direct transcription parameters
% pre-computed values for tau,w,D available for N = 10,20,50,100
p.ntA = 10; p.ntB = 10; % number of node points in Phase A and B

p.tauA = double(Gauss_nodes(p.ntA-1)); % scaled time horizon, Phase A
p.tauB = double(Gauss_nodes(p.ntB-1)); % scaled time horizon, Phase B
p.DA = double(Gauss_Dmatrix(p.tauA)); % differential approximation matrix (A)
p.DB = double(Gauss_Dmatrix(p.tauB)); % differential approximation matrix (B)
p.wA = double(Gauss_weights(p.tauA)); % for gaussian quadrature (A)
p.wB = double(Gauss_weights(p.tauB)); % for gaussian quadrature (B)
disp('Finished computing Gauss nodes and weights.')
% discretized variable indices in x = [y1A,y2A,uA,y1B,y2B,uB]:
nA = p.ntA; nB = p.ntB;
p.y1Ai = 1:nA; p.y2Ai = nA+1:2*nA; p.uAi = 2*nA+1:3*nA-1; m = 3*nA-1;
p.y1Bi = 1:nB; p.y2Bi = nB+1:2*nB; p.uBi = 2*nB+1:3*nB-1;
p.y1Bi = p.y1Bi + m; p.y2Bi = p.y2Bi + m; p.uBi = p.uBi + m;
p.tfi = 3*nB+m; % final time (free)
p.tsi = 3*nB+m+1; % switching time (free)
x0 = zeros(3*nB+m+1,1); % initial guess (all zeros)
x0(p.tfi) = 1;
x0(p.tsi) = 0.5;
x0(p.uAi) = -1;
x0(p.uBi) = 1;

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6,'Algorithm','sqp');

% solve the problem
[x,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);

% obtain the optimal solution
% Extract:
y1A = x(p.y1Ai); y2A = x(p.y2Ai); uA = x(p.uAi);
y1B = x(p.y1Bi); y2B = x(p.y2Bi); uB = x(p.uBi);
p.tf = x(p.tfi); % extract final time
p.ts = x(p.tsi); % extract switching time
% Unscale time horizon:
p.tA = (p.tauA*(p.ts-p.t0) + (p.ts+p.t0))/2;
p.tB = (p.tauB*(p.tf-p.ts) + (p.tf+p.ts))/2;
p.t = [p.tA; p.tB]; % Note that ts is not included here.

% obtain the costates
Lam1A = zeros(p.ntA+1, 1);
Lam2A = zeros(p.ntA+1, 1);
Lam1A(1) = lambda.eqnonlin(1);
Lam2A(1) = lambda.eqnonlin(2);
Lam1A(end) = -lambda.eqnonlin(3);
Lam2A(end) = -lambda.eqnonlin(5);
Lam1A(2:end-1) = lambda.eqnonlin(7:6+p.ntA-1)./p.wA + Lam1A(end);
Lam2A(2:end-1) = lambda.eqnonlin(6+p.ntA:6+2*p.ntA-2)./p.wA + Lam2A(end);
Lam1A = -Lam1A;
Lam2A = -Lam2A;
k = 6+2*p.ntA-2;
Lam1B = zeros(p.ntB+1, 1);
Lam2B = zeros(p.ntB+1, 1);
Lam1B(1) = -lambda.eqnonlin(3);
Lam2B(1) = -lambda.eqnonlin(5);
Lam1B(end) = -lambda.eqnonlin(4);
Lam2B(end) = -lambda.eqnonlin(6);
Lam1B(2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.wB + Lam1B(end);
Lam2B(2:end-1) = lambda.eqnonlin(k+p.ntB:k+2*p.ntB-2)./p.wB + Lam2B(end);
Lam1B = -Lam1B;
Lam2B = -Lam2B;
Lam.Lam1A = Lam1A; Lam.Lam2A = Lam2A;
Lam.Lam1B = Lam1B; Lam.Lam2B = Lam2B;

% calculate mu
muA = (2/(p.ts-p.t0))*lambda.ineqnonlin(1:p.ntA-1)./p.wA;
muB = (2/(p.tf-p.ts))*lambda.ineqnonlin(p.ntA:p.ntA+p.ntB-2)./p.wB;
Lam.muA = muA; Lam.muB = muB;

% plots
Plots_Gauss_BB2(y1A,y2A,uA,y1B,y2B,uB,Lam,p,'Pseudospectral')

% objective function
function J = objective(x,p)
    tf = x(p.tfi); % extract
    J = tf; % calculate objective
end

% constraint function
function [c,ceq] = constraints(x,p)
    % Extract:
    y1A = x(p.y1Ai); y2A = x(p.y2Ai); uA = x(p.uAi);
    y1B = x(p.y1Bi); y2B = x(p.y2Bi); uB = x(p.uBi);
    tf = x(p.tfi);
    ts = x(p.tsi);
    t0 = p.t0;
    
    % create matrices (p.nt x p.ns)
    YA = [y1A,y2A]; % For Gauss pseudospectral method, k=0 is needed here.
    YB = [y1B,y2B];
    FA = ((ts-t0)/2)*[y2A(2:end),uA]; % for Gauss PS, no k=0 here.
    FB = ((tf-ts)/2)*[y2B(2:end),uB];
    
    % initial state conditions
    ceq1 = y1A(1) - p.y10;
    ceq2 = y2A(1) - p.y20;

    ceq5A = p.DA*YA - FA; % differential equation approximation
    ceq5B = p.DB*YB - FB;
    
    c1A = abs(uA) - 1; % path constraints
    c1B = abs(uB) - 1;

    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq3B = y1B(1) + dot(p.wB,FB(:,1)) - p.y1f;
    ceq4B = y2B(1) + dot(p.wB,FB(:,2)) - p.y2f;
    
    % continuity requirement at switching time
    % necessary when implementing second phase
    ceq3A = y1A(1) + dot(p.wA,FA(:,1)) - y1B(1);
    ceq4A = y2A(1) + dot(p.wA,FA(:,2)) - y2B(1);
    
    c2 = ts - tf;
    c3A = max(uA) - min(uA) - 0.3;
    c3B = max(uB) - min(uB) - 0.3;
    
    % combine constraints
    c = [c1A;c1B;c2;c3A;c3B];
    ceq = [ceq1;ceq2;ceq3A;ceq3B;ceq4A;ceq4B;ceq5A(:);ceq5B(:)];
end


