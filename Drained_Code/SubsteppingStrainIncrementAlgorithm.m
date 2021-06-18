function [sigma1,PlasticWork1,DSTRAIN,flag,A,B] = SubsteppingStrainIncrementAlgorithm(dstrain1,PlasticWork0,m,b,testtype,sigma0,flag,A,B,MP)
Epois = MP(1); Em = MP(2); Elambda = MP(3); Tens = MP(4); Fm = MP(5); Feta = MP(6); Psi2 = MP(7); Pmu = MP(8);
Wc = MP(9); Wp = MP(10); Yh = MP(11); Yalpha = MP(12); Pa = MP(13); Psi1 = MP(14); Wrho = MP(15); Wd = MP(16);
ModelParameter = [Epois,Em,Elambda,Tens,Fm,Feta,Psi2,Pmu,Wc,Wp,Yh,Yalpha,Pa,Psi1,Wrho,Wd];
[f1,~,~,~,~,flag,A,B,~] = EvaluateFunctions(ModelParameter,sigma0,flag,A,B,PlasticWork0); % initial function value
% calculating the trial elastic stress or elastic shooting step
Ce = Elasticity(ModelParameter,sigma0); [DSTRAIN,dsigmat1] = EvaluateSigma(testtype,b,Ce,dstrain1); sigmat1 = sigma0+dsigmat1;
% evaluate the yield function for the trial elastic step
[f1t,f2,~,~,~,flag,A,B,~] = EvaluateFunctions(ModelParameter,sigmat1,flag,A,B,PlasticWork0); % after elastic shooting
% Check for elastic step
ftol = 5E-05;
if (f1t-f2) <= ftol
    sigmad1 = sigmat1; PlasticWorkd1 = PlasticWork0; % purely elastic step
    EDSTRAIN = DSTRAIN; PDSTRAIN = 0;
else
    if (f1-f2) <= ftol % the current step is elastic and plastic
        %elastic plus plastic step (determine alpha for the elastic segment)
        sigmat12 = sigma0; IntersectionTol = 5E-05;
        % algorithm to determine yield surface intersection parameter alpha
        for iter = 1:50
            sigmat11 = sigmat12;
            [f1t,~,~,~,~,flag,A,B,~] = EvaluateFunctions(ModelParameter,sigmat11,flag,A,B,PlasticWork0);
            [dfds,~] = DerivativeFunctions(MP,sigmat11); Ce = Elasticity(ModelParameter,sigma0);
            [DSTRAIN,dsigmat1] = EvaluateSigma(testtype,b,Ce,dstrain1);
            alpha = -(f1t-f2)./(dfds'*dsigmat1);
            sigmat12 = sigmat11+alpha*dsigmat1;
            if norm(sigmat12-sigmat11,2)/norm(sigmat11,2) < IntersectionTol
                break
            end
        end
        EDSTRAIN = alpha*DSTRAIN;
        dstrainplastic = (1-alpha)*dstrain1;
        sigmae1 = sigmat11;
    else
        % PURELY PLASTIC STEP
        %alpha = 0; purely plastic increment
        alpha = 0; %#ok<NASGU>
        dstrainplastic = dstrain1;
        sigmae1 = sigma0;
        EDSTRAIN = zeros(6,1);
    end
    
    %% forward Euler algorithm using m small steps to cover up dstrainplastic such that in each step consistency
    %  condition is satisfies
    
    incstrain = dstrainplastic/m; sigmad1 = sigmae1; plasticworkd1 = PlasticWork0;
    PDSTRAIN = zeros(6,1);
    for incdrift = 1:m
        % evaluate Cep and idlambda for given increment
        sigmad0 = sigmad1; plasticworkd0 = plasticworkd1;
        [~,~,~,~,~,flag,A,B,DFDWP0] = EvaluateFunctions(ModelParameter,sigmad0,flag,A,B,plasticworkd0);
        [Cep] = EvaluateCep(MP,sigmad0,DFDWP0);
        
        [dPDSTRAIN,idsigma] = EvaluateSigma(testtype,b,Cep,incstrain);
        [idlambda] = EvaluateLambda(MP,sigmad0,DFDWP0,dPDSTRAIN);
        [~,dgds] = DerivativeFunctions(MP,sigmad0);
        idplasticWork = idlambda*dgds'*sigmad0;
        PDSTRAIN = PDSTRAIN+dPDSTRAIN;
        sigmad1 = sigmad0+idsigma; plasticworkd1 = plasticworkd0+idplasticWork;
    end
    PlasticWorkd1 = plasticworkd1;
  
    %%
    % Check for the yield surface drift
    [f1d1,f2d1,~,~,~,flag,A,B] = EvaluateFunctions(ModelParameter,sigmad1,flag,A,B,plasticworkd1);
    sigmadt1 = sigmad1; plasticworkdt1 = plasticworkd1;
    if (f1d1-f2d1) > ftol
        for betaiter = 1:500
            sigmadt0 = sigmadt1; plasticworkdt0 = plasticworkdt1;
            [f1dt0,f2dt0,~,~,~,flag,A,B,dfdwt0] = EvaluateFunctions(ModelParameter,sigmadt0,flag,A,B,plasticworkdt0);
            if abs(f1dt0-f2dt0) < ftol
                break
            end
            Ce = Elasticity(MP,sigmadt0);
            [dfds,dgds] = DerivativeFunctions(MP,sigmadt0);
            beta = (f1dt0-f2dt0)/(dfds'*Ce*dgds-dfdwt0*dgds'*sigmadt0);
            plasticworkdt1 = plasticworkdt0+beta*dgds'*sigmadt0;
            sigmadt1 = sigmadt0-beta*Ce*dgds;
            betaiter
        end
        PlasticWorkd1 = plasticworkdt0;
        sigmad1 = sigmadt0;
    end
end
sigma1 = sigmad1;
PlasticWork1 = PlasticWorkd1;
DSTRAIN = EDSTRAIN+PDSTRAIN;
end
