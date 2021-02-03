function A = EMB_V2A(V, theta)
%EMB_V2A computes bilayer area from volumes
Rd=(V*3/4/pi).^(1/3);

Rl=Rd(:,1:end-1);
Rr=Rd(:,2:end);
C=cos(2*theta);
D = sqrt(Rl.^2+Rr.^2+2*Rl.*Rr.*C);
Rb = sqrt(-D.^4+2*D.^2.*Rl.^2+2.*D.^2.*Rr.^2-Rl.^4+2.*Rl.^2.*Rr.^2-Rr.^4)./(2*D);

A = pi*Rb.^2;
end

