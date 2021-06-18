function [I1,I2,I3,J2] = Invariant(sigma)
% Invariants of stress tensor 
% I1 = tr(sigma); I1 = sigma(1,1)+sigma(2,2)+sigma(3,3)
I1 = sigma(1,1)+sigma(2,1)+sigma(3,1);
% I2 = tr(cof(sigma)); I2 = 0.5*(I1^2-tr(sigma^2));
I2 = sigma(4,1)^2+sigma(5,1)^2+sigma(6,1)^2-(sigma(1,1)*sigma(2,1)+...
    sigma(2,1)*sigma(3,1)+sigma(3,1)*sigma(1,1));
% I3 = det(sigma); 
I3 = sigma(1,1)*(sigma(2,1)*sigma(3,1)-sigma(4,1)^2)-sigma(6,1)*...
    (sigma(6,1)*sigma(3,1)-sigma(5,1)*sigma(4,1))+sigma(5,1)*...
    (sigma(6,1)*sigma(4,1)-sigma(5,1)*sigma(2,1));
% J2 is the second invariant of the deviatoric tensor.
J2 = ((sigma(1,1)-sigma(2,1))^2+(sigma(2,1)-sigma(3,1))^2+...
    (sigma(3,1)-sigma(1,1))^2)/6+sigma(4,1)^2+sigma(5,1)^2+sigma(6,1)^2;