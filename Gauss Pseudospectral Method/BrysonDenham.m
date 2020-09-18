%--------------------------------------------------------------------------
% BrysonDenham.m
% Solve the Bryson-Denham problem using the differential Gauss
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
p.ns = 2; p.nu = 1; % number of states and controls
p.t0 = 0; p.tf = 1; % time horizon
p.y10 = 0; p.y1f = 0; p.y20 = 1; p.y2f = -1; % boundary conditions
p.l = 1/9;

% direct transcription parameters
% pre-computed values for tau,w,D available for N = 10,20,50,100
p.nt = 10; % number of node points, including -1
% p.nt = 20; % number of node points, including -1
% p.nt = 50; % number of node points, including -1
p.tau = double(Gauss_nodes(p.nt-1)); % scaled time horizon
p.D =  double(Gauss_Dmatrix(p.tau)); % differential approximation matrix
p.w = double(Gauss_weights(p.tau)); % for gaussian quadrature
% discretized variable indices in x = [y1,y2,u];
p.y1i = 1:p.nt; p.y2i = p.nt+1:2*p.nt; p.ui = 2*p.nt+1:3*p.nt-1;
x0 = zeros(p.nt*(p.ns+p.nu) - p.nu,1); % initial guess (all zeros)

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6,'Algorithm','sqp');

% solve the problem
x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);
% obtain the optimal solution
y1 = x(p.y1i); y2 = x(p.y2i); u = x(p.ui); % extract
p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon
% plots
Plots_Gauss_BD(y1,y2,u,p,'Pseudospectral')

% objective function
function J = objective(x,p)
    u = x(p.ui); % extract
    L = (u.^2)./2; % integrand
    J = ((p.tf-p.t0)/2)*dot(p.w,L); % calculate objective
end

% constraint function
function [c,ceq] = constraints(x,p)
    y1 = x(p.y1i); y2 = x(p.y2i); u = x(p.ui); % extract
    
    % create matrices (p.nt x p.ns)
    Y = [y1,y2]; % For Gauss pseudospectral method, k=0 is needed here.
    F = ((p.tf-p.t0)/2)*[y2(2:end),u]; % for Gauss PS, no k=0 here.
    
    % initial state conditions
    ceq1 = y1(1) - p.y10;
    ceq2 = y2(1) - p.y20;

    ceq5 = p.D*Y - F; % differential equation approximation
    c1 = y1 - p.l; % path constraints

    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq3 = y1(1) + dot(p.w,F(:,1)) - p.y1f;
    ceq4 = y2(1) + dot(p.w,F(:,2)) - p.y2f;
    
    % combine constraints
    c = c1; ceq = [ceq1;ceq2;ceq3;ceq4;ceq5(:);];
end