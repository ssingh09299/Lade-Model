function [DSTRAIN,dsigma] = EvaluateSigma(testtype,b,Cep,dstrain1)
if testtype == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for constant I1 test
    % Given the axial strain increment the other lateral strain increment
    % should be such that the 1) triaxial condition - dsigma2 = 0,
    % dsigma3 = 0; 2) constant I1 test - dsigma1 + dsigma2 + dsigma3 = 0; with
    % keeping b = constant which also means dsigma2 = b*dsigma1 + (1-b)*dsigma3
    
    SMAT1 = [1,1,1;b,-1,1-b]*[Cep(1,2),Cep(2,2),Cep(3,2);...
        Cep(1,3),Cep(2,3),Cep(3,3)]';
    SMAT2 = (-1)*[1,1,1;b,-1,1-b]*[Cep(1,1),Cep(2,1),Cep(3,1)]';
    
    % strain increment dstrain2 and dstrain3 is given by dstrain23
    dstrain(1,1) = dstrain1;  dstrain23 = (SMAT1\SMAT2)*dstrain1;
    dstrain(2,1) = dstrain23(1,1); dstrain(3,1) = dstrain23(2,1);
    dstrain(4,1) = 0.0; dstrain(5,1) = 0.0; dstrain(6,1) = 0.0;
    
    dSigma = Cep*dstrain;
elseif testtype == -1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for isotropic compression test
    % Given the axial strain increment the other lateral strain increment
    % should be such that the dsigma1 - dsigma2 = 0; dsigma1 - dsigma3 = 0;
    SMAT1 = [Cep(1,2)-Cep(2,2),Cep(1,3)-Cep(2,3);Cep(1,2)-Cep(3,2),Cep(1,3)-Cep(3,3)];
    SMAT2 = [Cep(1,1)-Cep(2,1);Cep(1,1)-Cep(3,1)];
    
    % strain increment dstrain2 and dstrain3 is given by dstrain23
    dstrain(1,1) = dstrain1;  dstrain23 = -(SMAT1\SMAT2)*dstrain1;
    dstrain(2,1) = dstrain23(1,1); dstrain(3,1) = dstrain23(2,1);
    dstrain(4,1) = 0.0; dstrain(5,1) = 0.0; dstrain(6,1) = 0.0;
    
     dSigma = Cep*dstrain;
elseif testtype == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for conventional triaxial compression or tension test
    % Given the axial strain increment the other lateral strain increment are
    % such that dsigma2 = 0; dsigma3 = 0;
    
    SMAT1 = [Cep(2,2),Cep(2,3);Cep(3,2),Cep(3,3)];
    SMAT2 = -[Cep(2,1);Cep(3,1)];
    
    % strain increment dstrain2 and dstrain3 is given by dstrain23
    dstrain(1,1) = dstrain1;  dstrain23 = (SMAT1\SMAT2)*dstrain1;
    dstrain(2,1) = dstrain23(1,1); dstrain(3,1) = dstrain23(2,1);
    dstrain(4,1) = 0.0; dstrain(5,1) = 0.0; dstrain(6,1) = 0.0;
    
    dSigma = Cep*dstrain;
    
else
    string = ['Error! The test type is neither the triaxial'...
        'nor the constant I1 test'];
    disp(string)
    string1 = ['Change it to either 1 or 0 and run it or other option would'...
        'be modify it'];
    error(string1)
end
dsigma = dSigma;
DSTRAIN = dstrain;
end
