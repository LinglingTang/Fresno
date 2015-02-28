function [ fval ] = solveFabts( M,L,U,P,Q,theU,theC,theV,t,s )
% solve M_inv(t,s)
n=size(M,1);
es=sparse(n,1);
es(s,1)=1;
[ M_inv_col_s ] = solveMue( L,U,P,Q,es );
val1=full(M_inv_col_s(t));
% the second part
tmp1=solveMue(L,U,P,Q,theU(:,1));
tmp2=solveMue(L,U,P,Q,theU(:,2));
B_invU=[tmp1,tmp2];
K=inv(theC)+theV'*(B_invU);
K_inv=inv(K);
tmp3=B_invU*K_inv;
tmp4=tmp3(t,:)*theV';
val2=-tmp4*M_inv_col_s;
fval=val1+val2;
end