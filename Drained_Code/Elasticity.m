function Ce = Elasticity(MP,sigma0)
% calculation for the parameters
Epois = MP(1); Em = MP(2); Elambda = MP(3); 
Pa = MP(13); 
% calculate the invariants
[I1,~,~,J2] = Invariant(sigma0);
% Modulus of elasticity given by Lade and Nelson (1989)
E = Em*Pa*((I1/Pa)^2+6*((1+Epois)/((1-2*Epois)))*(J2/Pa^2))^Elambda;
% Elastiicty tensor given by Lade and Jackobsen
Ce = [1-Epois,Epois,Epois,0,0,0;
    Epois,1-Epois,Epois,0,0,0;
    Epois,Epois,1-Epois,0,0,0;
    0,0,0,(1-2*Epois)/2,0,0;
    0,0,0,0,(1-2*Epois)/2,0;
    0,0,0,0,0,(1-2*Epois)/2]*(E/((1+Epois)*(1-2*Epois)));
