function [Ce,Cep,f1,f2,A,B,flag] = ElasticPlastic(sigma,Wc,Pa,Wp,Fm,Feta,Yalpha,Yh,Psi1,Wd,Wrho,Psi2,Pmu,Em,Epois,Elambda,flag,A,B)

W=Wc*(I1/Pa)^(Wp); 
fn = (I1^3/I3-27)*(I1/Pa)^(Fm); sl = fn/Feta; q = (Yalpha*sl)/(1-(1-Yalpha)*sl);
f1 = ((I1/Pa)^(Yh))*exp(q)*(Psi1*I1^3/I3-I1^2/I2);
f2=(W/(Wd*Pa))^(1/Wrho); 
g = ((Psi1*I1^3/I3-I1^2/I2+Psi2)*(I1/Pa)^(Pmu));
f = f1-f2;

if f < ftol
Cep = Ce;
else
if flag==0
    W=(f1^(Wrho))*(Pa*Wd);
    DFDWP = -(((1/(Pa*Wd))^(1/Wrho))*(1/Wrho))*(W)^(1/Wrho-1);
else
    W=-Pa*log(f1/A)/B;
    DFDWP = (A*B/Pa)*exp(-B*W/Pa);
end

if flag~=1 && sl<=1
    f2=(W/(Wd*Pa))^(1/Wrho);
elseif flag ~= 1 && sl>=1
    flag=1;
    B = Pa/(Wrho*W);
    A = f1*exp(B*W/Pa);
else
    f2=A*exp(-B*(W/Pa));
end

H=-Pmu*g*DFDWP;

DI1DSIGMA = [1,1,1,0,0,0]';
DI2DSIGMA = [-(sigma(2,1)+sigma(3,1)),-(sigma(3,1)+sigma(1,1)),-(sigma(1,1)+sigma(2,1)),2*sigma(4,1),2*sigma(5,1),2*sigma(6,1)]';
DI3DSIGMA = [(sigma(2,1)*sigma(3,1)-sigma(4,1)^2),(sigma(3,1)*sigma(1,1)-sigma(5,1)^2),(sigma(1,1)*sigma(2,1)-sigma(6,1)^2),...
    2*(sigma(6,1)*sigma(5,1)-sigma(1,1)*sigma(4,1)),2*(sigma(4,1)*sigma(6,1)-sigma(2,1)*sigma(5,1)),2*(sigma(5,1)*sigma(4,1)-sigma(3,1)*sigma(6,1))]';

DGDI1=(3+Pmu)*g/I1+(I1/I2-3*Psi2/I1)*(I1/Pa)^(Pmu); DGDI2=((I1/I2)^2)*(I1/Pa)^(Pmu); DGDI3=-Psi1*(I1^3/I3^2)*(I1/Pa)^(Pmu);
DQDI1=(Yalpha/(Feta*(1-(1-Yalpha)*sl)^2))*((Fm*sl*Feta/I1)+3*(I1^2/I3)*(I1/Pa)^Fm); DQDI3=-(Yalpha/(Feta*(1-(1-Yalpha)*sl)^2))*(I1^3/I3^2)*(I1/Pa)^Fm;
DFDI1=f1*(DQDI1+(3+Yh)/I1)+(I1/I2)*exp(q)*(I1/Pa)^Yh; DFDI2=((I1/I2)^2)*exp(q)*(I1/Pa)^(Yh); DFDI3=-Psi1*exp(q)*(I1^3/I3^2)*(I1/Pa)^(Yh)+f1*DQDI3;

DFDSIGMA=DFDI1*DI1DSIGMA+DFDI2*DI2DSIGMA+DFDI3*DI3DSIGMA; DGDSIGMA=DGDI1*DI1DSIGMA+DGDI2*DI2DSIGMA+DGDI3*DI3DSIGMA;
DFDS=DFDSIGMA; DGDS=DGDSIGMA;

% Elasto-plastic continuum tangent modulus
Cep = Ce-((Ce*DGDS)*(DFDS'*Ce))/(H+DFDS'*Ce*DGDS);
end