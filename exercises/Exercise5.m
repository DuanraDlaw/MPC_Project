%% Introduction:
clear all, close all;
set(0,'defaulttextinterpreter','latex');

% Parameters
dim = 2;
A = [0.7115 -0.4345; 0.4345 0.8853];
B = [0.2173; 0.0573];
Bd = zeros(dim,1);
C = [0 1];
Cd = 1;
d_r = 0.2;
M = [1; -1];
m = [3; 3];

%% Exercise 1: Observer design

% Compute the augmented state matrix
A_ = [A Bd; zeros(1,dim) 1];
B_ = [B; 0];
C_ = [C Cd];
p = [0.5, 0.6, 0.7];
K = place(A_', -C_', p);
L = K';

% Show that the estimate converges to the true values
x = [1; 2];
x_0 = [3; 0];
d_0 = 0;
dx = x - x_0;
dd = d_r - d_0;
d = zeros(3,dim);
d(:,1) = [dx; dd];
for i=1:30,
    d(:, i+1)=(A_ + L*C_)*d(:,i);
end

% Plots
figure(1)
hold on;
set(gca,'FontSize',14)
h1 = plot(d(1,:), d(2,:), 'Linestyle', '-', 'Marker', '+', 'color', 'b', 'LineWidth', 1);
h2 = plot(0, 0, 'Marker', 'o', 'color', 'g', 'LineWidth', 2);
xlabel('$\Delta x_1$', 'Fontsize', 18);
ylabel('$\Delta x_2$', 'Fontsize', 18);
hold off;

figure(2)
hold on;
set(gca,'FontSize',14)
h1 = stairs([0:30], d(3,:), 'Linestyle', '-', 'color', 'b', 'LineWidth', 1);
h2 = stairs([0:30], zeros(1,31), 'Linestyle', '--', 'color', 'g', 'LineWidth', 1);
xlabel('$t$', 'Fontsize', 18);
ylabel('$\Delta d$', 'Fontsize', 18);
hold off;

%% Exercise 2: Steady-state target computation

% Define optimization variables
x_s = sdpvar(2,1,'full');
u_s = sdpvar(1,1,'full');
d = sdpvar(1,1,'full');
r = sdpvar(1,1,'full');

con = [];
obj = 0;
con = [eye(dim)-A -B; C 0]*[x_s; u_s] == [Bd*d; r-Cd*d];
obj = u_s'*u_s;
con = [con, M*u_s <= m];

% Compile the matrices
SS_ctrl = optimizer(con, obj, [], [r d], [x_s; u_s]);
state = SS_ctrl{[1.2, 0.2]};

%% Exercise 3: MPC tracking

% Additional parameters
N = 5; % horizon
Nsim = 50; % number of simulation steps
Q = eye(dim); % stage cost (state)
R = 1; % stage cost (input)

% Compute terminal weight
Qf = dlyap(A,Q);

% Define optimization variables
dx = sdpvar(2,N,'full');
du = sdpvar(1,N,'full');
u_s = sdpvar(1, N, 'full');

% Define constraints and objective
con = [];
obj = 0;
for i = 1:N-1,
    con = [con, dx(:,i+1) == A*dx(:,i) + B*du(:,i)];      % System dynamics
    con = [con, M*du(:,i) <= m - M*u_s(:,i)];                        % Input constraints
    obj = obj + dx(:,i)'*Q*dx(:,i) + du(:,i)'*R*du(:,i);   % Cost function
end
obj = obj + dx(:,N)'*Qf*dx(:,N);  % Terminal constraint
    
% Compile the matrices
MPC_ctrl = optimizer(con, obj, [], [dx(:,1); u_s(:,1)], du(:,1));

% Initialize variables
x_r = zeros(2,Nsim+1);
x = zeros(2,Nsim+1);
d = zeros(1,Nsim+1);
y = zeros(1,Nsim+1);
u = zeros(1,Nsim);
x_r(:,1) = [1; 2];
x(:,1) = [3; 0]; % initial estimate
d(:,1) = 0; % initial estimate
y(:,1) = C*x_r(:,1) + Cd*d_r; % real y0
% r = 1;
r = 0.5;


% Simulation
for i = 1:Nsim,    
    % Obtain steady-state target
    SState = SS_ctrl{[r, d(:,i)]};
    x_s = SState(1:2);
    u_s = SState(3);
    
    % Initial state
    dx = x(:,i) - x_s;

    % Solve MPC problem to find the command u
    du = MPC_ctrl{[dx; u_s]};
    u(:,i) = du + u_s;
    
    % Evolution of real system
    x_r(:,i+1) = A*x_r(:,i) + B*u(:,i);
    y(:,i+1) = C*x_r(:,i+1) + Cd*d_r;
    
    % Estimate state and disturbance
    estimator = A_*[x(:,i);d(:,i)] + B_*u(:,i) + L*(C_*[x(:,i); d(:,i)] - y(:,i));
    x(:,i+1) = estimator(1:2);
    d(:,i+1) = estimator(3);
end

%% Plots

% The estimates converge to the true values
set(0,'defaulttextinterpreter','latex');

figure(4)
hold on;
set(gca,'FontSize',14)
h1 = plot(x(1,:), x(2,:), 'Linestyle', '-', 'Marker', '+', 'color', 'b', 'LineWidth', 1);
xlabel('$\hat{x_1}$', 'Fontsize', 18);
ylabel('$\hat{x_2}$', 'Fontsize', 18);
h2 = plot(x_r(1,:), x_r(2,:), 'Linestyle', '--', 'Marker', '+', 'color', 'g', 'LineWidth', 1);
xlabel('$x_1$', 'Fontsize', 18);
ylabel('$x_2$', 'Fontsize', 18);
legend([h1 h2], {'$\hat{x}$' , '$x$'}, 'Interpreter', 'latex', 'Fontsize', 18);
hold off;

figure(5)
hold on;
set(gca,'FontSize',14)
h1 = stairs([0:Nsim], d(1,:), 'Linestyle', '-', 'color', 'b', 'LineWidth', 1);
h2 = stairs([0:Nsim], d_r.*ones(1,Nsim+1), 'Linestyle', '--', 'color', 'g', 'LineWidth', 1);
xlabel('$t$', 'Fontsize', 18);
ylabel('$d$', 'Fontsize', 18);
legend([h1 h2], {'$\hat{d}$' , '$d$'}, 'Interpreter', 'latex', 'Fontsize', 18);
hold off;

% The output converges to the reference
figure(6)
hold on;
set(gca,'FontSize',14)
h1 = stairs([0:Nsim], y(1,:), 'Linestyle', '-', 'color', 'b', 'LineWidth', 1);
h2 = stairs([0:Nsim], r.*ones(1,Nsim+1), 'Linestyle', '--', 'color', 'g', 'LineWidth', 1);
xlabel('$t$', 'Fontsize', 18);
ylabel('$y$', 'Fontsize', 18);
legend([h1 h2], {'$y$' , '$r$'}, 'Interpreter', 'latex', 'Fontsize', 18);
hold off;

% The input does not violate the constraints
figure(7);
h1 = stairs([0:Nsim-1], u(1,:), 'LineWidth', 1);
set(gca,'FontSize',14)
set(gca,'XTick',[0:5:Nsim-1]);
set(gca,'YTick',[-4:1:4]);
xlim([0;Nsim-1]);
hold on;
h2 = plot([0:Nsim-1], 3*ones(1,Nsim),'Linestyle', '--', 'LineWidth', 1);
h3 = plot([0:Nsim-1], -3*ones(1,Nsim),'Linestyle', '--', 'LineWidth', 1);
xlabel('$t$', 'Fontsize', 18);
ylabel('$u$', 'Fontsize', 18);
legend([h1 h2 h3], {'$u$' , '$u_{high}$', '$u_{low}$'}, 'Interpreter', 'latex', 'Fontsize', 18);
hold off;