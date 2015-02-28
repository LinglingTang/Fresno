% CS 4220
% Porject 1
% Hang Chu(hc772) and Lingling Tang(lt375) @cornell.edu

% Setting up the data
load('roadNet-CA.mat');
tic
alpha = 0.9;
A = Problem.A;
n = length(A);
d = full(sum(A))';
I = find(d);
n = length(I);
d = d(I);
A = A(I,I);
D = spdiags(d, 0, n, n);
Dinv = spdiags(1./d, 0, n, n);
T = A * Dinv;
N = D - alpha * A;
M = speye(n) - alpha * T;

% Task 1
I = speye(n);
e1 = I(:,1);
tic;
u = M \ e1;
timeval = toc;
disp(['Backslash takes ',num2str(timeval), ' seconds']);

% Task 2
tic;
[L,U,P,Q] = lu(M);
timeval = toc;
disp(['LU decomposition takes ',num2str(timeval), ' seconds']);

% Task 4
tic;
[ newu ] = solveMue( L,U,P,Q,e1 );
timeval = toc;
disp(['Solving Mu=e_1 with LU decomposition takes ',num2str(timeval), ' seconds']);
diff = sum(abs(newu-u));
disp(['L1 distance between backslash solution and our solution is ',num2str(diff)]);

% Task 5
[ Mupdate,valid_update ] = computeMupdate( M,alpha,A,d,1,2 );

% Task 6
tic;
[ theU,theC,theV ] = svds( Mupdate );
timeval = toc;
disp(['Sparse SVD takes ',num2str(timeval), ' seconds']);
theU = theU(:,1:2);
theC = [theC(1,1),0;0,theC(2,2);];
theV =theV(:,1:2);
tic;
[ fval ] = solveFabts( M,L,U,P,Q,theU,theC,theV,180,1 );
timeval = toc;
disp(['Solving f(a,b)_ts takes ',num2str(timeval), ' seconds']);