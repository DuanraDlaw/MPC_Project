%%
%--------------------------------------------------------------------------
%                           Compute Invariant Sets
%--------------------------------------------------------------------------
clear all;
close all;

alpha = pi/6;
beta = 0.8;
A = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)].*beta;
H = [cos(pi/3) sin(pi/3); -cos(pi/3) -sin(pi/3); sin(pi/3) -cos(pi/3); -sin(pi/3) cos(pi/3)];
h = [2 1 2 5]';


X = Polyhedron(H,h);

% Plot the maximum invariant set OInf
S = X;

while true,
    F = S.A;
    f = S.b;
    preS = Polyhedron(F*A, f);
    prevS = S;
    S = Polyhedron([preS.A; prevS.A] , [preS.b; prevS.b]);
    if S == prevS,
        break
    end
end
OInf = S;
figure(1);
plot(X,'color', 'orange', OInf);
xlabel('x_1');
ylabel('x_2');
legend('X', 'O_{\infty}');

% Plot a trajectory where x0 is not in OInf and there exists an xi not in X

N = 30; % Size of the trajectory
x = zeros(2,10);
x(:,1) = [-3; 3];
for i=1:N,
    x(:,i+1)=A*x(:,i);
end
figure(2);
plot(X, 'color', 'orange', OInf);
xlabel('x_1');
ylabel('x_2');
hold on;
plot(x(1,:), x(2,:), 'Linewidth', 1.5);
legend('X', 'O_{\infty}', 'T_X');

% Plot several (>5) trajectories starting from various states within OInf,
% demonstrating that the entire trajectory {xi} remains within OInf

x1(:,1) = [-3 1];
x2(:,1) = [-2 2];
x3(:,1) = [-1 0];
x4(:,1) = [-2 0];
x5(:,1) = [0 2];
x6(:,1) = [2 1];

for i=1:N,
    x1(:,i+1)=A*x1(:,i);
    x2(:,i+1)=A*x2(:,i);
    x3(:,i+1)=A*x3(:,i);
    x4(:,i+1)=A*x4(:,i);
    x5(:,i+1)=A*x5(:,i);
    x6(:,i+1)=A*x6(:,i);
end
figure(3);
plot(X, 'color', 'orange', OInf);
xlabel('x_1');
ylabel('x_2');
legend('X', 'O_{\infty}');
hold on
plot(x1(1,:), x1(2,:), 'Linewidth', 1.5);
plot(x2(1,:), x2(2,:), 'Linewidth', 1.5);
plot(x3(1,:), x3(2,:), 'Linewidth', 1.5);
plot(x4(1,:), x4(2,:), 'Linewidth', 1.5);
plot(x5(1,:), x5(2,:), 'Linewidth', 1.5);
plot(x6(1,:), x6(2,:), 'Linewidth', 1.5);

%%
%--------------------------------------------------------------------------
%                   Compute Controlled Invariant Sets
%--------------------------------------------------------------------------


alpha = pi/6;
beta = 0.8;
A = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)].*beta;
B = [0.5 0.5]';
G = [1 -1]'; % slide 59 (constraint on u)
g = [0.5 0.5]'; % slide 59 (constraint on u)
H = [cos(pi/3) sin(pi/3); -cos(pi/3) -sin(pi/3); sin(pi/3) -cos(pi/3); -sin(pi/3) cos(pi/3)];
h = [2 1 2 5]';

X = Polyhedron(H,h);

% Plot the maximum controlled invariant set CInf
S = X;

while true,
    F = S.A;
    f = S.b;
    preS = projection(Polyhedron([F*A F*B; zeros(2,2) G], [f;g]), 1:2);
    prevS = S;
    S = Polyhedron([preS.A; prevS.A] , [preS.b; prevS.b]);
    if S == prevS,
        break
    end
end
CInf = S;
figure(1);
plot(X, 'color', 'orange', 'Alpha', 1, CInf, 'color', 'g', 'Alpha', 0.8);
hold on;


% Compute the optimal LQR controller K for Q = I, R = 1
Q = eye(2,2);
R = 1;

[P, ~, K] = dare(A,B,Q,R); % Use the infinite-horizon LQR optimization
K = -K;

% Plot the maximum invariant set OInf_ for the closed-loop system (u = Kx)
X_ = Polyhedron([H;G*K], [h;g]);
S = X_;

while true,
    F = S.A;
    f = S.b;
    preS = Polyhedron(F*(A+B*K), f);
    prevS = S;
    S = Polyhedron([preS.A; prevS.A] , [preS.b; prevS.b]);
    if S == prevS,
        break
    end
end
OInf_ = S;
plot(X_, 'color', 'c', 'Alpha', 0.8, OInf_, 'color', 'b', 'Alpha', 0.8);
xlabel('x_1');
ylabel('x_2');
legend('X', 'C_{\infty}', 'X_C', 'O_{\infty}');




