%--------------------------------------------------------------------------
% Rotating2.m
% (Using two-phase approach)
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
% qf = (1/sqrt(2))*[3/sqrt(14) -2/sqrt(14) 1/sqrt(14) 1]; % Quaternion final condition
w0 = [0 0 0]; % Angular velocity initial condition
wf = [0 0 0]; % Angular velocity final condition
% State boundary conditions, four quaternions, three angular velocity:
p.q10 = q0(1); p.q20 = q0(2); p.q30 = q0(3); p.q40 = q0(4);
p.q1f = qf(1); p.q2f = qf(2); p.q3f = qf(3); p.q4f = qf(4);
p.w10 = w0(1); p.w20 = w0(2); p.w30 = w0(3);
p.w1f = wf(1); p.w2f = wf(2); p.w3f = wf(3);

% direct transcription parameters
% pre-computed values for tau,w,D available for N = 10,20,50,100
p.ntA = 20; p.ntB = 20; % number of node points in Phase A and B

p.tauA = double(Gauss_nodes(p.ntA-1)); % scaled time horizon, Phase A
p.tauB = double(Gauss_nodes(p.ntB-1)); % scaled time horizon, Phase B
p.DA = double(Gauss_Dmatrix(p.tauA)); % differential approximation matrix (A)
p.DB = double(Gauss_Dmatrix(p.tauB)); % differential approximation matrix (B)
p.WA = double(Gauss_weights(p.tauA)); % for gaussian quadrature (A)
p.WB = double(Gauss_weights(p.tauB)); % for gaussian quadrature (B)
disp('Finished computing Gauss nodes and weights.')
% discretized variable indices: x = [q1A,q2A,...,w1A,...,u1A,...,q1B...];
nA = p.ntA; nB = p.ntB;
p.q1Ai = 1:nA; p.q2Ai = nA+1:2*nA;
p.q3Ai = 2*nA+1:3*nA; p.q4Ai = 3*nA+1:4*nA; m = 4*nA;
p.w1Ai = 1:nA; p.w2Ai = nA+1:2*nA; p.w3Ai = 2*nA+1:3*nA;
p.w1Ai = p.w1Ai + m; p.w2Ai = p.w2Ai + m;
p.w3Ai = p.w3Ai + m; m = m + 3*nA;
p.u1Ai = 1:nA-1; p.u2Ai = nA:2*nA-2; p.u3Ai = 2*nA-1:3*nA-3;
p.u1Ai = p.u1Ai + m; p.u2Ai = p.u2Ai + m;
p.u3Ai = p.u3Ai + m; m = m + 3*nA - 3;
p.q1Bi = 1:nB; p.q2Bi = nB+1:2*nB;
p.q3Bi = 2*nB+1:3*nB; p.q4Bi = 3*nB+1:4*nB;
p.q1Bi = p.q1Bi + m; p.q2Bi = p.q2Bi + m;
p.q3Bi = p.q3Bi + m; p.q4Bi = p.q4Bi + m; m = m + 4*nB;
p.w1Bi = 1:nB; p.w2Bi = nB+1:2*nB; p.w3Bi = 2*nB+1:3*nB;
p.w1Bi = p.w1Bi + m; p.w2Bi = p.w2Bi + m;
p.w3Bi = p.w3Bi + m; m = m + 3*nB;
p.u1Bi = 1:nB-1; p.u2Bi = nB:2*nB-2; p.u3Bi = 2*nB-1:3*nB-3;
p.u1Bi = p.u1Bi + m; p.u2Bi = p.u2Bi + m;
p.u3Bi = p.u3Bi + m; m = m + 3*nB - 3;
p.tfi = m+1; % final time (free)
p.tsi = m+2; % switching time (free)

% generate a reasonable initial guess
[qA, wA, uA, qB, wB, uB, tf, ts] = initialGuess(p, q0', qf');
x0 = zeros(m + 2,1);
x0(p.q1Ai) = qA(1,:); x0(p.q2Ai) = qA(2,:);
x0(p.q3Ai) = qA(3,:); x0(p.q4Ai) = qA(4,:);
x0(p.w1Ai) = wA(1,:); x0(p.w2Ai) = wA(2,:); x0(p.w3Ai) = wA(3,:);
x0(p.u1Ai) = uA(1,:); x0(p.u2Ai) = uA(2,:); x0(p.u3Ai) = uA(3,:);
x0(p.q1Bi) = qB(1,:); x0(p.q2Bi) = qB(2,:);
x0(p.q3Bi) = qB(3,:); x0(p.q4Bi) = qB(4,:);
x0(p.w1Bi) = wB(1,:); x0(p.w2Bi) = wB(2,:); x0(p.w3Bi) = wB(3,:);
x0(p.u1Bi) = uB(1,:); x0(p.u2Bi) = uB(2,:); x0(p.u3Bi) = uB(3,:);
x0(p.tfi) = tf; x0(p.tsi) = ts;

% Optimization algorithm options
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e6,'Algorithm','sqp');


% solve the problem; obtain the optimal solution
[x,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);
% x = x0;

% Extract:
q1A = x(p.q1Ai); q2A = x(p.q2Ai); q3A = x(p.q3Ai); q4A = x(p.q4Ai);
w1A = x(p.w1Ai); w2A = x(p.w2Ai); w3A = x(p.w3Ai);
u1A = x(p.u1Ai); u2A = x(p.u2Ai); u3A = x(p.u3Ai);
qA = [q1A, q2A, q3A, q4A];
wA = [w1A, w2A, w3A];
uA = [u1A, u2A, u3A];
q1B = x(p.q1Bi); q2B = x(p.q2Bi); q3B = x(p.q3Bi); q4B = x(p.q4Bi);
w1B = x(p.w1Bi); w2B = x(p.w2Bi); w3B = x(p.w3Bi);
u1B = x(p.u1Bi); u2B = x(p.u2Bi); u3B = x(p.u3Bi);
qB = [q1B, q2B, q3B, q4B];
wB = [w1B, w2B, w3B];
uB = [u1B, u2B, u3B];
p.tf = x(p.tfi); % extract final time
p.ts = x(p.tsi); % extract switching time
% Unscale time horizon:
p.tA = (p.tauA*(p.ts-p.t0) + (p.ts+p.t0))/2;
p.tB = (p.tauB*(p.tf-p.ts) + (p.tf+p.ts))/2;
p.t = [p.tA; p.tB]; % Note that ts is not included here.

% obtain the costates
LamqA = zeros(4, p.ntA+1);
LamwA = zeros(3, p.ntA+1);
LamqA(:,1) = lambda.eqnonlin(1:4);
LamwA(:,1) = lambda.eqnonlin(5:7);
LamqA(:,end) = -lambda.eqnonlin(8:11);
LamwA(:,end) = -lambda.eqnonlin(12:14);
k = 21;
LamqA(1,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamqA(1,end); k=k+p.ntA-1;
LamqA(2,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamqA(2,end); k=k+p.ntA-1;
LamqA(3,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamqA(3,end); k=k+p.ntA-1;
LamqA(4,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamqA(4,end); k=k+p.ntA-1;
LamwA(1,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamwA(1,end); k=k+p.ntA-1;
LamwA(2,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamwA(2,end); k=k+p.ntA-1;
LamwA(3,2:end-1) = lambda.eqnonlin(k+1:k+p.ntA-1)./p.WA + LamwA(3,end); k=k+p.ntA-1;
LamqA = -LamqA;
LamwA = -LamwA;
LamqB = zeros(4, p.ntB+1);
LamwB = zeros(3, p.ntB+1);
LamqB(:,1) = -lambda.eqnonlin(8:11);
LamwB(:,1) = -lambda.eqnonlin(12:14);
LamqB(:,end) = -lambda.eqnonlin(15:18);
LamwB(:,end) = -lambda.eqnonlin(19:21);
LamqB(1,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamqB(1,end); k=k+p.ntB-1;
LamqB(2,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamqB(2,end); k=k+p.ntB-1;
LamqB(3,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamqB(3,end); k=k+p.ntB-1;
LamqB(4,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamqB(4,end); k=k+p.ntB-1;
LamwB(1,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamwB(1,end); k=k+p.ntB-1;
LamwB(2,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamwB(2,end); k=k+p.ntB-1;
LamwB(3,2:end-1) = lambda.eqnonlin(k+1:k+p.ntB-1)./p.WB + LamwB(3,end); k=k+p.ntB-1;
LamqB = -LamqB;
LamwB = -LamwB;
Lam.LamqA = LamqA; Lam.LamwA = LamwA;
Lam.LamqB = LamqB; Lam.LamwB = LamwB;

% calculate mu
muA = zeros(3, p.ntA-1);
muB = zeros(3, p.ntB-1); b = 0;
muA(1,:) = (2/(p.ts-p.t0))*lambda.ineqnonlin(b+1:b+p.ntA-1)./p.WA; b = b+p.ntA-1;
muB(1,:) = (2/(p.tf-p.ts))*lambda.ineqnonlin(b+1:b+p.ntB-1)./p.WB; b = b+p.ntB-1;
muA(2,:) = (2/(p.ts-p.t0))*lambda.ineqnonlin(b+1:b+p.ntA-1)./p.WA; b = b+p.ntA-1;
muB(2,:) = (2/(p.tf-p.ts))*lambda.ineqnonlin(b+1:b+p.ntB-1)./p.WB; b = b+p.ntB-1;
muA(3,:) = (2/(p.ts-p.t0))*lambda.ineqnonlin(b+1:b+p.ntA-1)./p.WA; b = b+p.ntA-1;
muB(3,:) = (2/(p.tf-p.ts))*lambda.ineqnonlin(b+1:b+p.ntB-1)./p.WB; b = b+p.ntB-1;
Lam.muA = muA; Lam.muB = muB;

% plots
Plots_Gauss_Rotation2(qA,wA,uA,qB,wB,uB,Lam,p,'Pseudospectral')

% objective function
function J = objective(x,p)
    % Extract:
    tf = x(p.tfi);
    
    % Calculate objective:
    J = tf; % for minimum time
end

% constraint function
function [c,ceq] = constraints(x,p)
    % Extract:
    q1A = x(p.q1Ai); q2A = x(p.q2Ai); q3A = x(p.q3Ai); q4A = x(p.q4Ai);
    w1A = x(p.w1Ai); w2A = x(p.w2Ai); w3A = x(p.w3Ai);
    u1A = x(p.u1Ai); u2A = x(p.u2Ai); u3A = x(p.u3Ai);
    qA = [q1A, q2A, q3A, q4A];
    wA = [w1A(2:end), w2A(2:end), w3A(2:end)];
    uA = [u1A, u2A, u3A];
    q1B = x(p.q1Bi); q2B = x(p.q2Bi); q3B = x(p.q3Bi); q4B = x(p.q4Bi);
    w1B = x(p.w1Bi); w2B = x(p.w2Bi); w3B = x(p.w3Bi);
    u1B = x(p.u1Bi); u2B = x(p.u2Bi); u3B = x(p.u3Bi);
    qB = [q1B, q2B, q3B, q4B];
    wB = [w1B(2:end), w2B(2:end), w3B(2:end)];
    uB = [u1B, u2B, u3B];
    tf = x(p.tfi); % extract final time
    ts = x(p.tsi); % extract switching time
    t0 = p.t0;
    
    zerA = zeros(length(w1A),1);
    zerB = zeros(length(w1B),1);
    
    % Omega matrix
    Om1A = [zerA  w3A -w2A  w1A];
    Om2A = [-w3A zerA  w1A  w2A];
    Om3A = [ w2A -w1A zerA  w3A];
    Om4A = [-w1A -w2A -w3A zerA];
    Om1B = [zerB  w3B -w2B  w1B];
    Om2B = [-w3B zerB  w1B  w2B];
    Om3B = [ w2B -w1B zerB  w3B];
    Om4B = [-w1B -w2B -w3B zerB];
    
    % Dynamic equations
    dq1A = diag(0.5*Om1A*qA');
    dq2A = diag(0.5*Om2A*qA');
    dq3A = diag(0.5*Om3A*qA');
    dq4A = diag(0.5*Om4A*qA');
    dq1B = diag(0.5*Om1B*qB');
    dq2B = diag(0.5*Om2B*qB');
    dq3B = diag(0.5*Om3B*qB');
    dq4B = diag(0.5*Om4B*qB');
    
    dwA = p.jInv*(uA' - cross(wA',p.j*wA'));
    dwA = dwA';
    dwB = p.jInv*(uB' - cross(wB',p.j*wB'));
    dwB = dwB';
    
    % create matrices: (p.nt x p.ns) for Y, one fewer row for F.
    % For Gauss pseudospectral, k=0 is needed here:
    YA = [q1A,q2A,q3A,q4A,w1A,w2A,w3A];
    YB = [q1B,q2B,q3B,q4B,w1B,w2B,w3B];
    % For Gauss PS, no k=0 here:
    F0A = [dq1A(2:end),dq2A(2:end),dq3A(2:end),dq4A(2:end),dwA];
    F0B = [dq1B(2:end),dq2B(2:end),dq3B(2:end),dq4B(2:end),dwB];
    FA = ((ts-t0)/2)*F0A;
    FB = ((tf-ts)/2)*F0B;
    
    % initial state conditions
    ceq1 = q1A(1) - p.q10;
    ceq2 = q2A(1) - p.q20;
    ceq3 = q3A(1) - p.q30;
    ceq4 = q4A(1) - p.q40;
    ceq5 = w1A(1) - p.w10;
    ceq6 = w2A(1) - p.w20;
    ceq7 = w3A(1) - p.w30;
    
    % final state conditions
    % necessary for Gauss pseudospectral method
    ceq8B  = q1B(1) + dot(p.WB,FB(:,1)) - p.q1f;
    ceq9B  = q2B(1) + dot(p.WB,FB(:,2)) - p.q2f;
    ceq10B = q3B(1) + dot(p.WB,FB(:,3)) - p.q3f;
    ceq11B = q4B(1) + dot(p.WB,FB(:,4)) - p.q4f;
    ceq12B = w1B(1) + dot(p.WB,FB(:,5)) - p.w1f;
    ceq13B = w2B(1) + dot(p.WB,FB(:,6)) - p.w2f;
    ceq14B = w3B(1) + dot(p.WB,FB(:,7)) - p.w3f;
    
    % continuity requirement at switching time
    % necessary when implementing second phase
    ceq8A  = q1A(1) + dot(p.WA,FA(:,1)) - q1B(1);
    ceq9A  = q2A(1) + dot(p.WA,FA(:,2)) - q2B(1);
    ceq10A = q3A(1) + dot(p.WA,FA(:,3)) - q3B(1);
    ceq11A = q4A(1) + dot(p.WA,FA(:,4)) - q4B(1);
    ceq12A = w1A(1) + dot(p.WA,FA(:,5)) - w1B(1);
    ceq13A = w2A(1) + dot(p.WA,FA(:,6)) - w2B(1);
    ceq14A = w3A(1) + dot(p.WA,FA(:,7)) - w3B(1);

    ceq15A = p.DA*YA - FA; % differential equation approximation
    ceq15B = p.DB*YB - FB;
    ceq16A = sqrt(q1A.^2 + q2A.^2 + q3A.^2 + q4A.^2) - 1; % quaternion condition
    ceq16B = sqrt(q1B.^2 + q2B.^2 + q3B.^2 + q4B.^2) - 1;
    
    % path constraints
    c1A = abs(u1A) - p.u1max;
    c2A = abs(u2A) - p.u2max;
    c3A = abs(u3A) - p.u3max;
    c1B = abs(u1B) - p.u1max;
    c2B = abs(u2B) - p.u2max;
    c3B = abs(u3B) - p.u3max;
    c4 = -ts;
    c5 = ts - tf;
    c6A = max(u1A) - min(u1A) - 0.3*p.u1max;
    c7A = max(u2A) - min(u2A) - 0.3*p.u2max;
    c8A = max(u3A) - min(u3A) - 0.3*p.u3max;
    c6B = max(u1B) - min(u1B) - 0.3*p.u1max;
    c7B = max(u2B) - min(u2B) - 0.3*p.u2max;
    c8B = max(u3B) - min(u3B) - 0.3*p.u3max;
    
    % combine constraints
    c = [c1A;c1B;c2A;c2B;c3A;c3B;c4;c5; ...
%         c6A;c6B;c7A;c7B;c8A;c8B ...
        ];
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7; ...
        ceq8A;ceq9A;ceq10A;ceq11A;ceq12A;ceq13A;ceq14A; ...
        ceq8B;ceq9B;ceq10B;ceq11B;ceq12B;ceq13B;ceq14B; ...
        ceq15A(:);ceq15B(:);ceq16A;ceq16B];
end


% Find the initial guess.
function [qA, wA, uA, qB, wB, uB, tf, ts] = initialGuess(p, q0, qf)
    % Assumes the initial and final angular velocities are zero.
    % Also, q0 and qf must be column vectors.
    
    % number of time increments:
    M = 10000;
    
    % max control torque:
    Nmax = min([p.u1max, p.u2max, p.u3max]);
    
    q = q0; % current quaternion, varied
    qRem = quatMult(qf, quatInv(q)); % compute remaining rotation quaternion
    
    % First, get an order-of-magnitude estimate on the time for rotation:
    phiOri = 2*acos(qRem(4)); % angle in radians
    jApprox = max(eig(p.j));
    alpha = Nmax/jApprox;
    T = 2*sqrt(phiOri/alpha); disp(T)
    
    DeltaT = 1.2*T/M; % increment of time
    qArr = zeros(4,M); % array of quaternion components
    qArr(:,1) = q0;
    wArr = zeros(3,M); % array of angular velocity components
    w = wArr(:,1);
    NArr = zeros(3,M); % array of torque components
    t = 0; % time, varied
    
    % Compute initial axis of rotation:
    den = sqrt(qf(1)^2 + qf(2)^2 + qf(3)^2);
    un0 = qf(1:3)/den;
    un0Quat = qf; un0Quat(4) = 0; un0Quat(1:3) = un0;
    
    % temporary
    M = 1*M;
%     temp = zeros(M,1);
    
    % Compute Phase 1:
    for it = 1:M
        % Compute current axis of rotation:
        unQuat = un0Quat;
%         unQuat = quatMult(quatMult(quatConj(q),un0Quat),q);
%         unQuat = quatMult(quatMult(q,un0Quat),quatConj(q));
        un = unQuat(1:3);
        
        % Compute control torque:
        N = Nmax*un;
        NArr(:,it) = N;
        
        % Compute Omega matrix:
        wx = w(1); wy = w(2); wz = w(3);
        Om = [ 0   wz -wy  wx; ...
              -wz  0   wx  wy; ...
               wy -wx  0   wz; ...
              -wx -wy -wz  0 ];
        
        % Compute new quaternion value:
        DeltaQ = 0.5*DeltaT*Om*q;
        q = q + DeltaQ;
        qArr(:,it+1) = q;
        
        % Compute new angular velocity value:
        DeltaW = DeltaT*p.jInv*(N - cross(w,p.j*w));
        w = w + DeltaW;
        wArr(:,it+1) = w;
        
        t = t + DeltaT; % increment the time
        
        % Compute remaining rotation quaternion:
        qRem = quatMult(qf, quatInv(q));
        
        phiRem = 2*acos(qRem(4)); % angle of rotation remaining
%         temp(it) = phiRem*180/pi;
        
        if phiRem < phiOri/2
            itSwitch = it+1;
            break
        end
        
        % temporary
        if it == M
            itSwitch = M-1;
            disp('Something is wrong...')
        end
    end
    
%     plot(1:M,temp)
    
    % Extract quantities from Phase 1:
    disp('Initial guesses for ts and tf:')
    ts = t; disp(ts)
    tf = 2*ts; disp(tf)
    nA = p.ntA;
    qA = zeros(4, nA); wA = zeros(3, nA); uA = zeros(3, nA-1);
    m = (itSwitch - 1)/(nA - 1);
    for it = 1:nA-1
        currentIndex = (it - 1)*m + 1;
        currentIndex = round(currentIndex);
        qA(:,it) = qArr(:,currentIndex);
        wA(:,it) = wArr(:,currentIndex);
        uA(:,it) = NArr(:,currentIndex);
    end
    currentIndex = (currentIndex + itSwitch)/2;
    currentIndex = round(currentIndex);
    qA(:,nA) = qArr(:,currentIndex);
    wA(:,nA) = wArr(:,currentIndex);
    
    
    % Compute Phase 2:
    it = itSwitch;
    while ((t < tf) && (it < M))
        % Compute current axis of rotation:
        unQuat = un0Quat;
%         unQuat = quatMult(quatMult(quatConj(q),un0Quat),q);
%         unQuat = quatMult(quatMult(q,un0Quat),quatConj(q));
        un = unQuat(1:3);
        
        % Compute control torque:
        N = -Nmax*un; % note the negative sign for Phase 2
        NArr(:,it) = N;
        
        % Compute Omega matrix:
        wx = w(1); wy = w(2); wz = w(3);
        Om = [ 0   wz -wy  wx; ...
              -wz  0   wx  wy; ...
               wy -wx  0   wz; ...
              -wx -wy -wz  0 ];
        
        % Compute new quaternion value:
        DeltaQ = 0.5*DeltaT*Om*q;
        q = q + DeltaQ;
        qArr(:,it+1) = q;
        
        % Compute new angular velocity value:
        DeltaW = DeltaT*p.jInv*(N - cross(w,p.j*w));
        w = w + DeltaW;
        wArr(:,it+1) = w;
        
        t = t + DeltaT;
        it = it + 1;
    end
    
    % Extract quantities from Phase 2:
    itFinal = it; itBeg = itSwitch;
    nB = p.ntB;
    qB = zeros(4, nB); wB = zeros(3, nB); uB = zeros(3, nB-1);
    m = (itFinal - itBeg)/(nB - 1);
    for it2 = 1:nB-1
        currentIndex = (it2 - 1)*m + itBeg;
        currentIndex = round(currentIndex);
        qB(:,it2) = qArr(:,currentIndex);
        wB(:,it2) = wArr(:,currentIndex);
        uB(:,it2) = NArr(:,currentIndex);
    end
    currentIndex = (currentIndex + itFinal)/2;
    currentIndex = round(currentIndex);
    qB(:,nB) = qArr(:,currentIndex);
    wB(:,nB) = wArr(:,currentIndex);
    
end

