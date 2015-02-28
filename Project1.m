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
a = 1;
b = 2;
[ Mupdate,valid_update ] = computeMupdate( M,alpha,A,d,a,b );

% Task 6
t = 180;
s = 1;
tic;
[ theU,theC,theV ] = prepareUCV( n,a,b,A,alpha,d );
timeval = toc;
disp(['Decompositing UCV takes ',num2str(timeval), ' seconds']);
tic;
[ fval ] = solveFabts( M,L,U,P,Q,theU,theC,theV,t,s );
timeval = toc;
disp(['Solving f(a,b)_ts takes ',num2str(timeval), ' seconds']);

% Task 8
% find all roads that has degree higher than 1
allroad_sub = find(A);
[ allroad_ind1,allroad_ind2 ] = ind2sub( [n,n],allroad_sub );
degree1_points = find(d==1);
tf = ~ismember( allroad_ind2,degree1_points );
tmpind = find(tf);
allroad_ind2 = allroad_ind2(tmpind);
allroad_ind1 = allroad_ind1(tmpind);
tf = ~ismember( allroad_ind1,degree1_points );
tmpind = find(tf);
allroad_ind2 = allroad_ind2(tmpind);
allroad_ind1 = allroad_ind1(tmpind);
t = 180;
s = 1;
[ newL,newU,newP,newQ ] = lu(N);
et = sparse(n,1);
et(t) = 1;
[ w ] = solveMue( newL,newU,newP,newQ,et );
clear newL newU newP newQ
goodroad_ind1 = [];
goodroad_ind2 = [];
% find all roads that has upper bound higher than 1e-3
for ii=1:length(allroad_ind1)
    ii
    if ((allroad_ind1(ii)==t)&&(allroad_ind2(ii)==s)) || ((allroad_ind1(ii)==s)&&(allroad_ind2(ii)==t))
        continue;
    else
        [ upperbound ] = computeUpperBound( d,alpha,w,allroad_ind1(ii),allroad_ind2(ii),t );
        if upperbound > 1e-3
            goodroad_ind1 = [goodroad_ind1;allroad_ind1(ii)];
            goodroad_ind2 = [goodroad_ind2;allroad_ind2(ii)];
        end
    end
end
% find the optimal edge
tic;
es = sparse(n,1);
es(s) = 1;
[ us ] = solveMue( L,U,P,Q,es );
fval_base = full(us(t));
for ii=1:length(goodroad_ind1)
    [ upperbound ] = computeUpperBound( d,alpha,w,goodroad_ind1(ii),goodroad_ind2(ii),t );
    goodroad_upperbound(ii) = upperbound;
end
needchecking = ones(length(goodroad_ind1),1);
bestfval = -Inf;
besta = 0;
bestb = 0;
nowchecking = 0;
while sum(needchecking>0)
    for ii=1:length(goodroad_ind1)
        if needchecking(ii)==1
            nowchecking = ii;
            break;
        end
    end
    nowa = goodroad_ind1(nowchecking);
    nowb = goodroad_ind2(nowchecking);
    disp(['Checking ',int2str(nowa),' to ',int2str(nowb),'...']);
    [ theU,theC,theV ] = prepareUCV( n,nowa,nowb,A,alpha,d );
    [ fval ] = solveFabts( M,L,U,P,Q,theU,theC,theV,t,s );
    needchecking(nowchecking) = 0;
    if fval>bestfval
        bestfval = fval;
        besta = nowa; bestb = nowb;
    end
    for ii=1:length(goodroad_ind1)
        if (goodroad_ind2(ii)==nowa) && (goodroad_ind1(ii)==nowb)
            needchecking(ii) = 0;
        end
        if goodroad_upperbound(ii)<(bestfval-fval_base)
            needchecking(ii)=0;
        end
    end
end
timeval = toc;
disp(['Finding the optimal edge takes ',num2str(timeval), ' seconds']);
disp(['The optimal edge is ',int2str(besta),' to ',int2str(bestb),', with f(a,b) = ',num2str(bestfval)]);