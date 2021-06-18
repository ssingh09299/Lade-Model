function [dfds,dgds] = DerivativeFunctions(MP,sigma)
% Model parameters (MP)
Fm = MP(5); Feta = MP(6); Psi2 = MP(7); Pmu = MP(8);
Yh = MP(11); Yalpha = MP(12); Pa = MP(13); Psi1 = MP(14);

% Calculate the invariants
[I1,I2,I3,~] = Invariant(sigma);

% Derivatives of the invariants
% dI1/dsigma
DI1DSIGMA = [1,1,1,0,0,0]';
% dI2/dsigma
DI2DSIGMA = [-(sigma(2,1)+sigma(3,1)),-(sigma(3,1)+sigma(1,1)),...
    -(sigma(1,1)+sigma(2,1)),2*sigma(4,1),2*sigma(5,1),2*sigma(6,1)]';
% dI3/dsigma
DI3DSIGMA = [(sigma(2,1)*sigma(3,1)-sigma(4,1)^2),(sigma(3,1)*...
    sigma(1,1)-sigma(5,1)^2),(sigma(1,1)*sigma(2,1)-sigma(6,1)^2),...
    2*(sigma(6,1)*sigma(5,1)-sigma(1,1)*sigma(4,1)),2*(sigma(4,1)...
    *sigma(6,1)-sigma(2,1)*sigma(5,1)),2*(sigma(5,1)*sigma(4,1)-...
    sigma(3,1)*sigma(6,1))]';

% evalute the functons
fn = (I1^3/I3-27)*(I1/Pa)^(Fm);
% calculate sl,q,f1,g
sl = fn/Feta; q = (Yalpha*sl)/(1-(1-Yalpha)*sl);
f1 = ((I1/Pa)^(Yh))*exp(q)*(Psi1*I1^3/I3-I1^2/I2);
g = ((Psi1*I1^3/I3-I1^2/I2+Psi2)*(I1/Pa)^(Pmu));

% calculate the derivatives of plastic potential function
DGDI1 = (3+Pmu)*g/I1+(I1/I2-3*Psi2/I1)*(I1/Pa)^(Pmu);
DGDI2 = ((I1/I2)^2)*(I1/Pa)^(Pmu); DGDI3=-Psi1*(I1^3/I3^2)*(I1/Pa)^(Pmu);
DQDI1 = (Yalpha/(Feta*(1-(1-Yalpha)*sl)^2))*((Fm*sl*Feta/I1)+3*(I1^2/I3)*...
    (I1/Pa)^Fm);
% calculate the derivatives of the yield function
DQDI3 = -(Yalpha/(Feta*(1-(1-Yalpha)*sl)^2))*(I1^3/I3^2)*(I1/Pa)^Fm;
DFDI1 = f1*(DQDI1+(3+Yh)/I1)+(I1/I2)*exp(q)*(I1/Pa)^Yh;
DFDI2 = ((I1/I2)^2)*exp(q)*(I1/Pa)^(Yh);
DFDI3 = -Psi1*exp(q)*(I1^3/I3^2)*(I1/Pa)^(Yh)+f1*DQDI3;

DFDSIGMA = DFDI1*DI1DSIGMA+DFDI2*DI2DSIGMA+DFDI3*DI3DSIGMA;
DGDSIGMA = DGDI1*DI1DSIGMA+DGDI2*DI2DSIGMA+DGDI3*DI3DSIGMA;
dfds = DFDSIGMA;
dgds = DGDSIGMA;
