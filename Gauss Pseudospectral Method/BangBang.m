%--------------------------------------------------------------------------
% BangBang.m
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
% p.nt = 11; % number of node points, including -1
p.nt = 21; % number of node points, including -1
% p.nt = 51; % number of node points; this takes a few minutes
p.tau = double(Gauss_nodes(p.nt-1)); % scaled time horizon
p.D =  double(Gauss_Dmatrix(p.tau)); % differential approximation matrix
p.w = double(Gauss_weights(p.tau)); % for gaussian quadrature
% discretized variable indices in x = [y1,y2,u];
p.y1i = 1:p.nt; p.y2i = p.nt+1:2*p.nt; p.ui = 2*p.nt+1:3*p.nt-1;
p.tfi = 3*p.nt; % Be sure to add 1 in the zeros below:
x0 = zeros(p.nt*(p.ns+p.nu) - p.nu + 1,1); % initial guess (all zeros)
x0(p.nt*(p.ns+p.nu) - p.nu + 1) = 1;

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6,'Algorithm','sqp');

% solve the problem
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);
% obtain the optimal solution
y1 = x(p.y1i); y2 = x(p.y2i); u = x(p.ui); % extract
p.tf = x(p.tfi); % extract final time
p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon
% plots
Plots_Gauss_BB(y1,y2,u,p,'Pseudospectral')

% objective function
function J = objective(x,p)
    tf = x(p.tfi); % extract
    J = tf; % calculate objective
end

% constraint function
function [c,ceq] = constraints(x,p)
    y1 = x(p.y1i); y2 = x(p.y2i); u = x(p.ui); % extract
    tf = x(p.tfi);  % extract
    
    % create matrices (p.nt x p.ns)
    Y = [y1,y2]; % For Gauss pseudospectral method, k=0 is needed here.
    F = ((tf-p.t0)/2)*[y2(2:end),u]; % for Gauss PS, no k=0 here.
    
    % initial state conditions
    ceq1 = y1(1) - p.y10;
    ceq2 = y2(1) - p.y20;

    ceq5 = p.D*Y - F; % differential equation approximation
    c1 = abs(u) - 1; % path constraints

    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq3 = y1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq4 = y2(1) + dot(p.w,F(:,2)) - p.y2f;
    
    % combine constraints
    c = c1; ceq = [ceq1;ceq2;ceq3;ceq4;ceq5(:);];
end