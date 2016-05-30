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
N = 20;

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
    constraints = constraints + (-zpMax <= x(1,i+1) <= zpMax);
    constraints = constraints + (-angleMax <= x(2:3,i+1) <= angleMax);
    constraints = constraints + (-angledMax <= x(5:6,i+1) <= angledMax);
    constraints = constraints + (-gammadMax <= x(7,i+1) <= gammadMax);
   
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


%%
pause

%% Nonlinear model simulation - no disturbance
fprintf('Running the FIRST NL model simulation...\n')
sim('simulation1.mdl') 

%%
pause

%% Disturbance estimation
%estimator


%% Offset free MPC
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')

pause

%% Final simulation
fprintf('Running the FINAL NL model simulation...\n')
sim('simulation2.mdl') 
pause
%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')






