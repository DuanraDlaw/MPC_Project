% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance
load('quadData.mat')
outerController = getOuterController(Ac, 'cplex');
disp('Data successfully loaded')



%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%

N = 10;
T = 2;

x0 = [-1 10*pi/180 -10*pi/180 120*pi/180 0 0 0]';

Q = eye(7);
R = eye(4);

zpMax = 1;                  % [m/s]
alphaMax = 10*pi/180;       % [rad]
angleMax = [alphaMax; alphadMax];
alphadMax = 15*pi/180;      % [rad/s]
angledMax = [alphadMax; alphadMax];
gammadMax = 60*pi/180;      % [rad/s]

x = sdpvar(7,N+1,'full');
u = sdpvar(4,N,'full');

[~,S,~] = dlqr(sys.A, sys.B, Q, R, zeros(7, 4));

constraints = [];
for i = 1:N
    constraints = constraints + (x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i));
end
constraints = constraints + (x(:,1) == x0);

for i = 1:N
    constraints = constraints + (-zpMax <= x(1,i+1) <= zpMax);
    constraints = constraints + (-angleMax <= x(2:3,i+1) <= angleMax);
    constraints = constraints + (-angledMax <= x(4:5,i+1) <= angledMax);
    constraints = constraints + (-gammadMax <= x(6,i+1) <= gammadMax);
    
    constraints = constraints + (0 <= u(:,i) <= 1);
end

objective = 0;
for i = 1:N
    objective = objective +  x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end 
objective = objective + x(:,N+1)'*S*x(:,N+1);

options = sdpsettings('solver','qpip');
innerController = optimizer(constraints, objective, options, x(:,1), u(:,1));
simQuad(sys, innerController, x0, T);

%%
pause

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






