% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance

clear;
close all;
clc;

load('quadData.mat')
outerController = getOuterController(Ac, 'gurobi');
disp('Data successfully loaded')

%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%

% Horizon steps
N = 10;

% Parameters for cost function 
Q = diag([1e4 5e5 5e5 1 1 1 1]);
R = eye(4);

% Constraints
zpMax = 1;                  % [m/s]
alphaMax = 10*pi/180;       % [rad]
angleMax = [alphaMax; alphaMax];
alphadMax = 15*pi/180;      % [rad/s]
angledMax = [alphadMax; alphadMax];
gammadMax = 60*pi/180;      % [rad/s]

% Steady-state input constraints
uMin = -us;
uMax = ones(4,1)-us;

% Terminal constraint using dlqr
[~,S,~] = dlqr(sys.A, sys.B, Q, R, zeros(7, 4));

% System definition
x = sdpvar(7,N,'full'); u = sdpvar(4,N,'full');
constraints = []; objective = 0;
for i = 1:N-1
    % System dynamics
    constraints = constraints + (x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i));
   
    % State constraints
    constraints = constraints + (-zpMax <= x(1,i) <= zpMax);
    constraints = constraints + (-angleMax <= x(2:3,i) <= angleMax);
    constraints = constraints + (-angledMax <= x(5:6,i) <= angledMax);
    constraints = constraints + (-gammadMax <= x(7,i) <= gammadMax);
   
    % Input Constraints
    constraints = constraints + (uMin <= u(:,i) <= uMax);
    
    % Cost function
    objective = objective +  x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
objective = objective + x(:,N)'*S*x(:,N);

% Simulate the model
T = 10;% Time horizon
x0 = [-1 10*pi/180 -10*pi/180 120*pi/180 0 0 0]'; % Initial state

options = sdpsettings('solver','gurobi');
innerController = optimizer(constraints, objective, options, x(:,1), u(:,1));
simQuad(sys, innerController, x0, T);

%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n')
% yalmip('clear');
close all;

% Define steady-state parameters
C = [eye(4) zeros(4,3)];

% Define optimization variables
x_r = sdpvar(7,1,'full');
u_r = sdpvar(4,1,'full');
r = sdpvar(4,1,'full');

constraints = []; objective = 0;
% Steady-state
constraints = [eye(7)-sys.A -sys.B; C zeros(4,4)]*[x_r; u_r] == [zeros(7,1); r];

% State constraints
constraints = constraints + (-zpMax <= x_r(1) <= zpMax);
constraints = constraints + (-angleMax <= x_r(2:3) <= angleMax);
constraints = constraints + (-angledMax <= x_r(5:6) <= angledMax);
constraints = constraints + (-gammadMax <= x_r(7) <= gammadMax);

% Input Constraints
constraints = constraints + (uMin <= u(:,i) <= uMax);

% Objective function
% objective = x_r'*Q*x_r + u_r'*R*u_r;
objective = u_r'*R*u_r;

% Compile the matrices
SS_ctrl = optimizer(constraints, objective, options, r, [x_r; u_r]);
r = [-1 10*pi/180 -10*pi/180 120*pi/180]';
state = SS_ctrl{r};

%% --------------------------------------------------------------------------
% MPC Tracking 

% Define optimization variables
x = sdpvar(7,1,'full');
r = sdpvar(4,1,'full');
dx = sdpvar(7,N,'full');
du = sdpvar(4,N,'full');
u_r = zeros(4,1); % u_r is close to 0 for arbitrary references (previous part)
x_r = [r; zeros(3,1)];

% Define constraints and objective
constraints = [];
objective = 0;

% Initial condition
constraints = constraints + (dx(:,1) == x - x_r);
for i = 1:N-1,
    % System dynamics
    constraints = constraints + (dx(:,i+1) == sys.A*dx(:,i) + sys.B*du(:,i));   
    % State constraints
    constraints = constraints + (-zpMax - x_r(1) <= dx(1,i) <= zpMax - x_r(1));                    
    constraints = constraints + (-angleMax - x_r(2:3) <= dx(2:3,i) <= angleMax - x_r(2:3));
    constraints = constraints + (-angledMax - x_r(5:6) <= dx(5:6,i) <= angledMax - x_r(5:6));
    constraints = constraints + (-gammadMax - x_r(7) <= dx(7,i) <= gammadMax - x_r(7));
    % Input constraints
    constraints = constraints + (uMin - u_r <= du(:,i) <= uMax - u_r);                       
    % Cost function
    objective = objective + dx(:,i)'*Q*dx(:,i) + du(:,i)'*R*du(:,i);            
end
% constraints = constraints + (dx(:,N) == 0);
objective = objective + dx(:,N)'*S*dx(:,N);     % Terminal constraint

innerController = optimizer(constraints, objective, options, [x(:,1); r(:,1)], du(:,1));

%% Tracking of constant reference signal
ref = 0.5* [-1 10*pi/180 -10*pi/180 120*pi/180]';
x0 = zeros(7,1);
simQuad( sys, innerController, x0, T, ref);

%% Tracking of slowly varying reference signal

T = 20;
ref = 0.5* [-1 10*pi/180 -10*pi/180 120*pi/180]';
coefs = linspace(0.5, 1, T/sys.Ts+1);
c_ref = repmat(ref,1, T/sys.Ts+1);
var_ref = zeros(4, T/sys.Ts+1);
for i=1:T/sys.Ts+1,
    var_ref(:,i) = c_ref(:,i)*coefs(i);
end

x0 = zeros(7,1);
simQuad( sys, innerController, x0, T, var_ref);


%%
pause

%% Nonlinear model simulation - no disturbance
fprintf('Running the FIRST NL model simulation...\n')
sim('simulation1.mdl') 

%%
pause

%% Disturbance estimation
%estimator
% Compute the augmented state matrix
A_ = [sys.A eye(7); zeros(7) eye(7)];
B_ = [sys.B; zeros(7,4)];
C_ = [eye(7) zeros(7)];
p = [0.5*ones(7,1); 0.4*ones(7,1)]; % Desired poles
K = place(A_', C_', p);
L = K';

filter.Af = A_ - L*C_;
filter.Bf = [B_ L];

%% Offset free MPC
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')
yalmip('clear')

% Define optimization variables
x_hat = sdpvar(7,1,'full');
r = sdpvar(4,1,'full');
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
d_hat = sdpvar(7,1,'full'); % Constant mean of disturbance
% u_r = zeros(4,1); % u_r is close to 0 for arbitrary references
% x_r = [r; zeros(3,1)];

% Define constraints and objective
constraints = [];
objective = 0;

% Target condition with disturbance
steady_state = [eye(7)-sys.A, -sys.B; C, zeros(4,4)]\[eye(7)*d_hat; r-zeros(4,7)*d_hat];
x_r = steady_state(1:7);
u_r = steady_state(8:11);
% Initial condition
constraints = constraints + (x(:,1) == x_hat);
for i = 1:N-1,
    % System dynamics with disturbance
    constraints = constraints + (x(:,i+1)== sys.A*x(:,i) + sys.B*u(:,i) + d_hat);   
    % State constraints
    constraints = constraints + (-zpMax <= x(1,i) <= zpMax);                    
    constraints = constraints + (-angleMax <= x(2:3,i) <= angleMax);
    constraints = constraints + (-angledMax <= x(5:6,i) <= angledMax);
    constraints = constraints + (-gammadMax <= x(7,i) <= gammadMax);
    % Input constraints
    constraints = constraints + (uMin <= u(:,i) <= uMax);
    % Cost function
    objective = objective + (x(:,i)-x_r)'*Q*(x(:,i)-x_r) + (u(:,i)-u_r)'*R*(u(:,i) - u_r);            
end
% constraints = constraints + (x(:,N) == x_r);
objective = objective + (x(:,N)-x_r)'*S*(x(:,N)-x_r);     % Terminal constraint

innerController = optimizer(constraints, objective, options, [x_hat(:,1); r(:,1); d_hat], u(:,1));

%% Step reference tracking with disturbance
close all;
T = 20;
ref = [-0.5 5*pi/180 -5*pi/180 90*pi/180]';

simQuad(sys, innerController, x0, T, ref, filter);

%% Slowly-varying tracking with disturbance

T = 20;
ref = 0.5* [-1 10*pi/180 -10*pi/180 120*pi/180]';
coefs = linspace(0.5, 1, T/sys.Ts+1);
c_ref = repmat(ref,1, T/sys.Ts+1);
var_ref = zeros(4, T/sys.Ts+1);
for i=1:T/sys.Ts+1,
    var_ref(:,i) = c_ref(:,i)*coefs(i);
end

x0 = zeros(7,1);
simQuad( sys, innerController, x0, T, var_ref, filter);


%%
pause

%% Final simulation
fprintf('Running the FINAL NL model simulation...\n')
sim('simulation2.mdl') 

%%
pause
%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')






