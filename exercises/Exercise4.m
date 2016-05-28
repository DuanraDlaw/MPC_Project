%%
%--------------------------------------------------------------------------
%                           Implement MPC
%--------------------------------------------------------------------------
close all

% Initialize parameters
A = [0.9752 1.4544; -0.0327 0.9315];
B = [0.0248; 0.0327];
F = [1 0; -1 0; 0 1; 0 -1];
f = [5; 5; 0.2; 0.2];
M = [1;-1];
m = [1.75;1.75];

% USING CODE FROM PREVIOUS WEEK
% Compute the optimal LQR controller K for Q = I, R = 1
Q = 1*eye(2);
R = 20;

[P, ~, K] = dare(A,B,Q,R); % Use the infinite-horizon LQR optimization
K = -K;
% P = weight matrix
% K = controller gain

% Plot the maximum invariant set OInf_ for the closed-loop system (u = Kx)
plot(Polyhedron(F, f))
X_ = Polyhedron([F;M*K], [f;m]);
S = X_;

while true,
    Floop = S.A;
    floop = S.b;
    preS = Polyhedron(Floop*(A+B*K), floop);
    prevS = S;
    S = Polyhedron([preS.A; prevS.A] , [preS.b; prevS.b]);
    if S == prevS,
        break
    end
end
OInf_ = S;
set(0,'defaulttextinterpreter','latex');
figure(1);
set(gca,'FontSize',14)
plot(X_, 'color', 'c', 'Alpha', 0.8, OInf_, 'color', 'b', 'Alpha', 0.8);
xlabel('$x_1$', 'Fontsize', 14);
ylabel('$x_2$', 'Fontsize', 14);
legend('X_c', 'O_{\infty}');

% USING MPT3 LQRSYSTEM

sys = LTISystem('A', A, 'B', B);
sys.x.max = [5 0.2]';
sys.x.min = [-5 -0.2]';
sys.u.max = 1.75;
sys.u.min = -1.75;

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

K_mpt = sys.LQRGain();
S_mpt = sys.LQRSet();
P_mpt = sys.LQRPenalty().weight;

OInf_mpt = S_mpt;
figure(2);
set(gca,'FontSize',14)
plot(X_, 'color', 'c', 'Alpha', 0.8, OInf_mpt, 'color', 'b', 'Alpha', 0.8);
xlabel('$x_1$', 'Fontsize', 14);
ylabel('$x_2$', 'Fontsize', 14);
legend('X_c', 'O_{\infty}^{mpt}');

% Compute matrices for solving MPC with optimization function
N = 10;
x_0 = [3 0]';
Ff = S.A;
ff = S.b;
Qf = P;
T = [kron(diag(ones(N-1,1),-1), -A) + eye(N*length(A)) kron(eye(N),-B)];
t = kron(eye(N,1), A)*x_0;
G = blkdiag(kron(eye(N-1), F), Ff, kron(eye(N), M));
g = [kron(ones(N-1,1), f); ff; kron(ones(N,1), m)];
H = blkdiag(kron(eye(N-1), Q), Qf, kron(eye(N), R));

% [zopt, fval, flag] = quadprog(H, zeros(3*N,1), G, g, T, t);

% Simulate closed-loop system
Nsim = 100; % Number of simulation steps
x = zeros(2, Nsim+1);
x(:,1) = x_0;
u = zeros(1, Nsim);
for i=1:Nsim,
    t = kron(eye(N,1), A)*x(:,i);
    [zopt, fval, flag] = quadprog(H, zeros(3*N,1), G, g, T, t);
    display(flag == 1);
    x(1,i+1) = zopt(1);
    x(2,i+1) = zopt(2);
    u(1,i) = zopt(2*N+1);
end

% Plots
figure(3);
set(gca,'FontSize',14)
h1 = plot(Polyhedron(F,f), 'color', 'r', 'Alpha', 0.5);
hold on;
h2 = plot(x(1,:), x(2,:), 'color', 'k', 'LineWidth', 1);
xlabel('$x_1$', 'Fontsize', 14);
ylabel('$x_2$', 'Fontsize', 14);
legend([h1 h2], {'X' , 'x'});
hold off;
figure(4);
set(gca,'FontSize',14)
h1 = stairs([0:Nsim-1], u(1,:), 'LineWidth', 1);
set(gca,'XTick',[0:2:Nsim-1]);
set(gca,'YTick',[0:2:Nsim-1]);
xlim([0;Nsim-1]);
hold on;
h2 = plot([0:Nsim-1], 1.75*ones(1,Nsim), 'LineWidth', 1);
h3 = plot([0:Nsim-1], -1.75*ones(1,Nsim), 'LineWidth', 1);
xlabel('$t$', 'Fontsize', 14);
ylabel('$u$', 'Fontsize', 14);
legend([h1 h2 h3], {'u' , 'u_{high}', 'u_{low}'});
hold off;
figure(5);
set(gca,'FontSize',14)
hold on;
[hAx, ~, ~]= plotyy([0:Nsim], x(1,:), [0:Nsim], x(2,:));
xlabel('$t$', 'FontSize', 14);
ylabel(hAx(1), '$x_1$', 'FontSize', 14);
ylabel(hAx(2), '$x_2$', 'FontSize', 14);
hold off;

%%
%--------------------------------------------------------------------------
%                       Implement MPC using YALMIP
%--------------------------------------------------------------------------

% Define optimization variables
x = sdpvar(2,N,'full');
u = sdpvar(1,N,'full');

% Define constraints and objective
con = [];
obj = 0;
for i = 1:N-1,
    con = [con, x(:,i+1) == A*x(:,i) + B*u(:,i)];      % System dynamics
    con = [con, F*x(:,i) <= f];                        % State constraints
    con = [con, M*u(:,i) <= m];                        % Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);   % Cost function
end
con = [con, Ff*x(:,N) <= ff];   % Terminal constraint
obj = obj + x(:,N)'*Qf*x(:,N);  % Terminal constraint
    
% Compile the matrices
ctrl = optimizer(con, obj, [], x(:,1), u(:,1));

% Can now compute the optimal control input using
x = zeros(2, Nsim+1);
x(:,1) = x_0;
u = zeros(1, Nsim);
for i = 1:Nsim,
    [uopt, isfeasible] = ctrl{x(:,i)};
    u(i)=uopt(1);
    x(:,i+1) = A*x(:,i) + B*u(i);
end

% Plots
figure(6);
set(gca,'FontSize',14)
h1 = plot(Polyhedron(F,f), 'color', 'r', 'Alpha', 0.5);
hold on;
h2 = plot(x(1,:), x(2,:), 'color', 'k', 'LineWidth', 1);
xlabel('$x_1$', 'Fontsize', 14);
ylabel('$x_2$', 'Fontsize', 14);
legend([h1 h2], {'X' , 'x^{yalmip}'});
hold off;
figure(7);
set(gca,'FontSize',14)
h1 = stairs([0:Nsim-1], u(1,:), 'LineWidth', 1);
set(gca,'XTick',[0:2:Nsim-1]);
set(gca,'YTick',[0:2:Nsim-1]);
xlim([0;Nsim-1]);
hold on;
h2 = plot([0:Nsim-1], 1.75*ones(1,Nsim), 'LineWidth', 1);
h3 = plot([0:Nsim-1], -1.75*ones(1,Nsim), 'LineWidth', 1);
xlabel('$t$', 'Fontsize', 14);
ylabel('$u$', 'Fontsize', 14);
legend([h1 h2 h3], {'u^{yalmip}' , 'u_{high}', 'u_{low}'});
hold off;
figure(8);
set(gca,'FontSize',14)
hold on;
[hAx, hLine1, hLine2]= plotyy([0:Nsim], x(1,:), [0:Nsim], x(2,:));
xlabel('$t$', 'FontSize', 14);
ylabel(hAx(1), '$x_1$', 'FontSize', 14);
ylabel(hAx(2), '$x_2$', 'FontSize', 14);
hold off;
