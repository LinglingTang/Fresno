function [ Mupdate,valid_update ] = computeMupdate( M,alpha,A,d,a,b )
% Intuitively, to shut down a road, the road must exist.
% Computationally, if for some city there was only one road connected to
% it, then we shouldn't shut it down to avoid divided by zero.
n=length(d);
if (d(a)<=1)||(d(b)<=1)
    disp('Oh no! You cannot shut down this road!');
    valid_update=0;
    Mupdate=sparse(n,n);
else
    valid_update=1;
    d(a)=d(a)-1;
    d(b)=d(b)-1;
    Dinv=spdiags(1./d,0,n,n);
    A(a,b)=0;
    A(b,a)=0;
    T=A*Dinv;
    Mnew=speye(n)-alpha*T;
    Mupdate=Mnew-M;
end
end