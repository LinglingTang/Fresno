function [ theU,theC,theV ] = prepareUCV( n,a,b,A,alpha,d )
theU = zeros(n,2);
roadsa = find(A(:,a)); roadsb = find(A(:,b));
theU(b,1) = alpha / d(a);
for ii=1:length(roadsa)
    if roadsa(ii)==b
        continue;
    else
        theU(roadsa(ii),1)=alpha / d(a) - alpha / (d(a) - 1);
    end
end
theU(a,2) = alpha / d(b);
for ii=1:length(roadsb)
    if roadsb(ii)==a
        continue;
    else
        theU(roadsb(ii),2)=alpha / d(b) - alpha / (d(b) - 1);
    end
end
theC = eye(2);
theV = zeros(n,2);
theV(a,1) = 1; theV(b,2) = 1;
end