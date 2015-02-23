% CS 4220
% Porject 1

% Task 1
load('roadNet-CA.mat');
tic
alpha = 0.9;
n = length(Problem.A);
d = full(sum(Problem.A))';
T = Problem.A*spdiags(1./d, 0, n, n); 
I = speye(n);
M = I-alpha*T;
e1=I(:,1);
u=M\e1;
toc

% Task 2
[L,U,P,Q] = lu(M);

% Task 3
