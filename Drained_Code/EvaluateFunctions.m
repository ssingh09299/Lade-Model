function [f1,f2,g,sl,q,flag,A,B,DFDWP] = EvaluateFunctions(MP,sigma,flag,A,B,W)
% The function evaluates f1,f2,g,sl,q,flag,A,B for the given step, the
% command is [f1,f2,g,sl,q,flag,A,B] = EvaluateFunctions(MP,sigma,...
% flag,A,B,W)
Fm = MP(5); Feta = MP(6); Psi2 = MP(7); Pmu = MP(8); Yh = MP(11); 
Yalpha = MP(12); Pa = MP(13); Psi1 = MP(14); Wrho = MP(15);
Wd = MP(16);
% calculate the invariants and failure function fn
[I1,I2,I3,~] = Invariant(sigma);
fn = (I1^3/I3-27)*(I1/Pa)^(Fm); 
% calculate sl,q,f1,g
sl = fn/Feta; q = (Yalpha*sl)/(1-(1-Yalpha)*sl);
f1 = ((I1/Pa)^(Yh))*exp(q)*(Psi1*I1^3/I3-I1^2/I2);
g = ((Psi1*I1^3/I3-I1^2/I2+Psi2)*(I1/Pa)^(Pmu));
% calculate hardening or softening function based on the regime using
% regulator sl and flag
% if flag is not 1 and sl < 1 implies hardening
if flag~=1 && sl<1
    f2=(W/(Wd*Pa))^(1/Wrho);
    % flag remains not equal to 1 till peak is reached where sl = 1
    DFDWP = -(((1/(Pa*Wd))^(1/Wrho))*(1/Wrho))*(W)^(1/Wrho-1);
elseif flag ~= 1 && sl>=1
    % change the flag to 1 once it reaches peak
    DFDWP = -(((1/(Pa*Wd))^(1/Wrho))*(1/Wrho))*(W)^(1/Wrho-1);
    flag=1;
    B = Pa/(Wrho*W);
    A = f1*exp(B*W/Pa);
    f2= A*exp(-B*(W/Pa));
else
    % once flag is equal to 1, we are in softening regime
    f2=A*exp(-B*(W/Pa));
    DFDWP = (A*B/Pa)*exp(-B*W/Pa);
end


