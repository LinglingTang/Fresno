function [ upperbound ] = computeUpperBound( d,alpha,w,a,b,t )
upperbound=(1/d(t))*(1+alpha)/(1-alpha)*(full(w(a))/(d(a)-1)+(full(w(b))/(d(b)-1)));
end