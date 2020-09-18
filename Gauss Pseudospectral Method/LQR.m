%--------------------------------------------------------------------------
% LQR.m
% Solve the LQR problem using the differential Gauss
% pseudospectral method
% (namely Gauss nodes and Gaussian quadrature)
% Method developed by David Benson, MIT PhD thesis
% LQR problem from page 126.
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
p.ns = 1; p.nu = 1; % number of states and controls
p.t0 = 0; p.tf = 5; % time horizon
p.y10 = 1; p.y1f = 0; % boundary conditions

% direct transcription parameters
% pre-computed values for tau,w,D available for N = 10,20,50,100
p.nt = 10; % number of node points, including -1
% p.nt = 20; % number of node points, including -1
% p.nt = 50; % number of node points, including -1
p.tau = double(Gauss_nodes(p.nt-1)); % scaled time horizon
p.D =  double(Gauss_Dmatrix(p.tau)); % differential approximation matrix
p.w = double(Gauss_weights(p.tau)); % for gaussian quadrature
% discretized variable indices in x = [y1,u];
p.y1i = 1:p.nt; p.ui = p.nt+1:2*p.nt-1;
x0 = zeros(p.nt*(p.ns+p.nu) - p.nu,1); % initial guess (all zeros)

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6,'Algorithm','sqp');

% solve the problem
[x,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);

% obtain the optimal solution
y1 = x(p.y1i); u = x(p.ui); % extract
p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon

% obtain the costate
Lam = zeros(p.nt + 1, 1);
Lam(1) = lambda.eqnonlin(1);
Lam(end) = -lambda.eqnonlin(2);
Lam(2:end-1) = lambda.eqnonlin(3:end)./p.w + Lam(end);
Lam = -Lam;

% plots
Plots_Gauss_LQR(y1,u,Lam,p,'Pseudospectral')

% objective function
function J = objective(x,p)
    % extract
    u = x(p.ui);
    y1 = x(p.y1i);
    
    L = (y1(2:end).^2 + u.^2)./2; % integrand
    J = ((p.tf-p.t0)/2)*dot(p.w,L); % calculate objective
end
% constraint function
function [c,ceq] = constraints(x,p)
    y1 = x(p.y1i); u = x(p.ui); % extract
    
    % create matrices (p.nt x p.ns)
    Y = y1; % For Gauss pseudospectral method, k=0 is needed here.
    F = ((p.tf-p.t0)/2)*(y1(2:end) + u); % for Gauss PS, no k=0 here.
    
    % initial state conditions
    ceq1 = y1(1) - p.y10;

    ceq5 = p.D*Y - F; % differential equation approximation

    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq3 = y1(1) + dot(p.w,F(:,1)) - p.y1f;
    
    % combine constraints
    c = -1; ceq = [ceq1;ceq3;ceq5(:);];
end