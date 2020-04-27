function Area=Bilayer_area(a,b,c)

l1 = a; 
l2 = b; 
theta = c; 

A=l1./2;
B=l2./2;
C=cos(2*theta);

D = sqrt(A^2+B^2+2*A*B*C);

r = sqrt(-D^4+2*D^2*A^2+2*D^2*B^2-A^4+2*A^2*B^2-B^4)/(2*D);

Area = pi*r.^2;

end

