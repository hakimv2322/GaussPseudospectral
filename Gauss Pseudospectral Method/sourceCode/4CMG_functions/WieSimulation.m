% This recreates the simulation of Wie et al. (2001)
% "Singularity Robust Steering Logic for Redundant Single-Gimbal Control
% Moment Gyros"
% Journal Guidance, Control, and Dynamics. Vol. 24, No. 5.

% Setup parameters
beta = 53.13; % skew angle, in degrees
h = 1000; % each CMG momentum, in SI units (constant)
deltaMax = 2; % max gimbal rate, in rad/s
slewMax = 10; % max slew rate, in deg/s
slewMax = slewMax*pi/180; % convert to rad/s
J = diag([21400,20100,5000]); % spacecraft tensor of inertia, in SI units
Jinv = inv(J);
totalTime = 12; % total time of simulation, in seconds
dt = 1/10000; % interval of time, in seconds

% Control parameters
omegaN = 3; % in rad/s
zeta = 0.9;
T = 1e1; % time constant of integral control, in seconds (should be 10)
k = omegaN^2 + 2*zeta*omegaN/T;
c = 2*zeta*omegaN + 1/T;

% Initial/final conditions
% q: spacecraft quaternion, w: spacecraft angular velocity (body frame),
% delta: gimbal angles, deltaDot: gimbal angle rates
q0 = [0; 0; 0; 1];
qf = [0.4; 0; 0; sqrt(0.84)];
delta0 = [0; 0; 0; 0]; % Case 1
% delta0 = [90; 0; -90; 0]; % Case 2 (elliptic singularity)
delta0 = delta0*pi/180; % convert to rad
w0 = [0; 0; 0];
deltaDot1D = [0; 0; 0; 0]; % first derivative of deltaDot

% Quaternion error matrix (for Wie's way)
r1 = [qf(4), qf(3), -qf(2), -qf(1)];
r2 = [-qf(3), qf(4), qf(1), -qf(2)];
r3 = [qf(2), -qf(1), qf(4), -qf(3)];
r4 = [qf(1), qf(2), qf(3), qf(4)];
QEmatrix = [r1; r2; r3; r4];

% Quantities over time
tArray = 0:dt:totalTime;
N = length(tArray); % number of columns for each quantity
qArray = zeros(4,N);
wArray = zeros(3,N);
deltaArray = zeros(4,N);
deltaDotArray = zeros(4,N);

% Initialize
qArray(:,1) = q0;
wArray(:,1) = w0;
deltaArray(:,1) = delta0;
intE = [0; 0; 0]; % integral of error quaternion vector
counter = 0; % for one-time reset of integral


% Calculate maneuver simulation
for it = 1:N-1
    
    t = tArray(it);
    q = qArray(:,it);
    w = wArray(:,it);
    delta = deltaArray(:,it);
    deltaDot = deltaDotArray(:,it);
    qv = q(1:3);
    A = Jacobian4(h, beta, delta*180/pi);
    hv = angVector4(h, beta, delta*180/pi);
    
    % Calculate error quaternion vector
%     qErr = quatMult(qf, quatInv(q)); % my way
    qErr = QEmatrix*q; % Wie's way
    err = qErr(1:3);
    
    % Calculate saturated torque (tau) and CMG torque (u)
    limitValues = Limiter(A, err, J, slewMax, deltaMax, k, c);
    errTerm = err + intE/T;
    satErr = [0; 0; 0];
    for it2 = 1:3
        satErr(it2) = sat(errTerm(it2), -limitValues(it2), limitValues(it2));
    end
    tau = -J*(2*k*satErr + c*w);
    u = -tau - cross(w, hv);
    
    % Calculate generalized robust pseudoinverse (Asharp)
    A = A/h; % to normalize A
    lam = 0.01*exp(-10*det(A*A'));
    eps1 = 0.01*sin(0.5*pi*t);
    eps2 = 0.01*sin(0.5*pi*t + pi/2);
    eps3 = 0.01*sin(0.5*pi*t + pi);
    r1 = [1, eps3, eps2];
    r2 = [eps3, 1, eps1];
    r3 = [eps2, eps1, 1];
    E = [r1; r2; r3];
    Asharp = A'*inv(A*A' + lam*E);
    Asharp = Asharp/h; % to un-normalize Asharp
    A = A*h; % to un-normalize A
    
    % Calculate commanded gimbal angle rates
    deltaDotc = Asharp*u;
    mu = max(abs(deltaDotc));
    if mu >= deltaMax
        deltaDotc = deltaDotc*deltaMax/mu;
    end
    
    % Calculate derivatives
    wDot = -Jinv*(cross(w, J*w+hv) + A*deltaDot);
    qvDot = 0.5*(-cross(w, qv) + q(4)*w);
    q4Dot = -0.5*dot(w, qv);
    qDot = [qvDot; q4Dot];
    
    % Gimbal rate dynamics
    % Don't include gimbal dynamics:
    deltaDotArray(:,it+1) = deltaDotc;
    % Include gimbal dynamics:
%     deltaDotArray(:,it+1) = deltaDot + dt*deltaDot1D;
%     deltaDot2D = (50^2)*(deltaDotc - deltaDot) - 2*0.7*50*deltaDot1D;
%     deltaDot1D = deltaDot1D + dt*deltaDot2D;
    
    % Update quantities
    qArray(:,it+1) = q + dt*qDot;
    wArray(:,it+1) = w + dt*wDot;
    deltaArray(:,it+1) = delta + dt*deltaDot;
    
    % Calculate integral of error (Trapezoid Rule)
    if it == 1
        intEalt = err*dt/2;
    elseif max(abs(err)) < 0.01 && counter == 0
        % one-time reset of the integral when err gets small
        intEalt = err*dt/2; % toggle for trapezoid rule
%         intE = [0; 0; 0]; % toggle for rectangle rule
        counter = 1;
    else
        intEalt = intEalt + err*dt;
    end
    intE = intEalt - err*dt/2; % toggle for trapezoid rule
%     intE = intE + err*dt; % toggle for rectangle rule
end

% Save quantities
Wie.dt = dt;
Wie.tArray = tArray;
Wie.qArray = qArray;
Wie.wArray = wArray;
Wie.deltaArray = deltaArray;
Wie.deltaDotArray = deltaDotArray;

% Convert gimbal angles into degrees
deltaArray = deltaArray*180/pi;

% Convert angular velocity components into deg/s
wArray = wArray*180/pi;

% Plot quaternion components
% plot(tArray, qArray(1,:))
% title('q1')
% figure
% plot(tArray, qArray(2,:))
% title('q2')
% figure
% plot(tArray, qArray(3,:))
% title('q3')
% figure
% plot(tArray, qArray(4,:))
% title('q4')

% Plot angular velocity components
% figure
% plot(tArray, wArray(1,:))
% title('w1'); ylabel('deg/s'); xlabel('seconds')
% figure
% plot(tArray, wArray(2,:))
% title('w2'); ylabel('deg/s'); xlabel('seconds')
% figure
% plot(tArray, wArray(3,:))
% title('w3'); ylabel('deg/s'); xlabel('seconds')

% Plot gimbal angles
% figure
% plot(tArray, deltaArray(1,:))
% title('delta1'); ylabel('degrees'); xlabel('seconds')
% figure
% plot(tArray, deltaArray(2,:))
% title('delta2'); ylabel('degrees'); xlabel('seconds')
% figure
% plot(tArray, deltaArray(3,:))
% title('delta3'); ylabel('degrees'); xlabel('seconds')
% figure
% plot(tArray, deltaArray(4,:))
% title('delta4'); ylabel('degrees'); xlabel('seconds')

% Plot gimbal angle rates
% figure
% plot(tArray, deltaDotArray(1,:))
% title('d/dt delta1'); ylabel('rad/s'); xlabel('seconds')
% figure
% plot(tArray, deltaDotArray(2,:))
% title('d/dt delta2'); ylabel('rad/s'); xlabel('seconds')
% figure
% plot(tArray, deltaDotArray(3,:))
% title('d/dt delta3'); ylabel('rad/s'); xlabel('seconds')
% figure
% plot(tArray, deltaDotArray(4,:))
% title('d/dt delta4'); ylabel('rad/s'); xlabel('seconds')


