function  [dlambda] = EvaluateLambda(MP,sigma,DFDWP,DSTRAIN)
% The function evaluates f1,f2,g,sl,q,flag,A,B for the given step, the
% command is [f1,f2,g,sl,q,flag,A,B] = EvaluateFunctions(MP,sigma,...
% flag,A,B,W)
Psi2 = MP(7); Pmu = MP(8); 
Pa = MP(13); Psi1 = MP(14);
% calculate the invariants and failure function fn
[I1,I2,I3,~] = Invariant(sigma); 
% calculate g
g = ((Psi1*I1^3/I3-I1^2/I2+Psi2)*(I1/Pa)^(Pmu));
[DFDS,DGDS] = DerivativeFunctions(MP,sigma);
H=-Pmu*g*DFDWP;
% Elasticity tensor
Ce = Elasticity(MP,sigma);
% Lambda as given by elasto-plastic theory
dlambda = (Ce*DFDS)'*DSTRAIN/(H+DFDS'*Ce*DGDS);
