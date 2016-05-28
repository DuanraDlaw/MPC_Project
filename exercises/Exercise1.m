clear all;
close all;

% Main parameters
N = 200;
A = [4/3 -2/3; 1 0];
B = [1 0]';
C = [-2/3 1];

% Exercise 1:

% Weights
Q = C'*C + 0.001*eye(2);
R = 0.001;
P_f = Q;

% % For finite horizon using least-squares optimization
% B_ = kron(eye(N),B);
% A_ = kron(diag(ones(N-1,1),-1), A) - eye(N*length(A));
% C_ = kron(eye(N,1), -A);
% Q_ = kron(eye(N), Q);
% R_ = kron(eye(N), R);
% 
% F = -inv(A_)*B_;
% G = A_\C_;
% 
% K = -inv(R_ + F'*Q_*F)*F'*Q_*G;

% For finite horizon using dynamic programming
    % start at step N and compute the minimum of the cost function
        H = Q;
    % iterate backward (DP iteration)
    K = zeros(N-1,2);
    for i = (N-1):(-1):0
        K(i+1,:) = -inv(R + B'*H*B)*B'*H*A;
        H = Q + K(i+1,:)'*R*K(i+1,:) + (A+B*K(i+1,:))'*H*(A+B*K(i+1,:));
    end

% Exercise 2:
M = 50; % number of simulation steps
x = [10 10]';
s = zeros(2,M);
s(:,1) = x;
v = 1;
p = zeros(2,N,2);
p(:,1,1) = x;
p(:,1,2) = x;

for i=1:M,
    u = K*x;
    x = A*x + B*u(1);
    s(:,i+1) = x;
    if i==1 || i==2 ,
        p(:,1,v)= s(:,i);
        for k=1:N,
            p(:,k+1,v) = A*p(:,k,v) + B*u(k);
        end
        v=v+1;
    end
end

figure(1);
plot(p(1,:,1),p(2,:,1),'-x', p(1,:,2),p(2,:,2), '-x');

figure(2);
plot(s(1,:),s(2,:),'-x');
axis([-5 15 -5 15]);

% stable for N >= 7

% Exercise 3: 
% For infinite-horizon LQR:
[P, ~, G] = dare(A,B,P_f,R);
K = -G;

M = 50; % number of simulation steps
x = [10 10]';
s = zeros(2,M);
s(:,1) = x;

for i=1:M,
    u = K*x;
    x = A*x + B*u(1);
    s(:,i+1) = x;
end

figure(3);
plot(s(1,:),s(2,:),'-x');
axis([-5 15 -5 15]);

% Try to plot absolute distance to (0,0)as a function of time step

