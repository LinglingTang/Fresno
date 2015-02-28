function [ u ] = solveMue( L,U,P,Q,e )
u=P*e;
u=L\u;
u=U\u;
u=Q*u;
end