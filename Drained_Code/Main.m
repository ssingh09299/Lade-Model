%{
This is the main file where you input the model parameters and test
conditions. The code uses forward Euler substepping algorithm for
integration of the model.

The code is valid for any arbitrary stress path, however, in this file only
consolidation, conventional triaxial compression, and constant mean
effective stress path (for different intermediate principal stress ratios)
are implemented.

If you require to change these for your customized stress path you have to
understand the code. Or send a request to Saurabh Singh for customization.

If you found a bug in the code or you would like to contribute, please
contact Saurabh (contact address below).

Please keep "Lade Model Integration: Git Hub" as your subject line, if you
contact me.

Contributions:
This code is "solely" written by Saurabh Singh when he was a PhD student at
the department of Civil Engineering at Indian Institute of Science,
Bangalore.
Author: Saurabh Singh
Contact: saurabhhbti08@gmail.com

Request:
If you use this code, please 

License: GPL v.3
%}

%{
Refer to Lade's papers (or my thesis) for more information on the model.

Usage:
The input parameters are present in Parameters.xsls, please change the
parameters there.

Additionally, you would require to change the control parameters for
elemental test conditions for which model integration is performed.

The results are in the .mat file
%}


set(0, 'DefaultAxesFontSize',20)
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontSize', 20)
set(0,'defaultTextFontName', 'Times New Roman')

%% CONTROL PARAMETERS

% Parameter for integration
% the yield surface drift
m = 80;

% condtions to pick up from for consolidation_type and shearing_type
% constant I1 test - 1; conventional triaxial - 0; isotropic compression test - -1;

% CONDITION FOR CONSOLIDATION STAGE
consolidation_type = -1; iConsolidationPressure = 300;

% Number of steps for consolidation
consolidation_steps = 90000; 

% Initial mean effective stress after consolidation
isigma = 10;

% Axial strain increment during consolidation
dstrain1 = 5E-05;

% CONDITION OF SHEARING STAGE
shearing_type = 0; 

% Intermediate principal stress ratio for shearing
b = 1.0;

% Number of steps for shearing
shearing_steps = 300000;

% Axial strain increment during shearing
dstrain2 = 5E-07;

%% No user input required here onwards %%
% However you are strongly welcome to go through the code.

%% Input parameters - imported from Parameters.xlsx
parameters = xlsread('Parameters.xlsx',1,'B1:B13');

Epois = parameters(1); Em = parameters(2); Elambda = parameters(3);
Tens = parameters(4)*101.3; Fm = parameters(6); Feta = parameters(5);
% Plastic potential, work hardening and yield parameters
Psi2 = parameters(8); Pmu = parameters(9); Wc = parameters(10);
Wp = parameters(11); Yh = parameters(13); Yalpha = parameters(12);
% dependent parameters
Pa = 101.3; Psi1 = 0.00155*Fm^(-1.27); Wrho = Wp/Yh; Wd = Wc/(27*Psi1+3)^(Wrho);

ModelParameter = [Epois,Em,Elambda,Tens,Fm,Feta,Psi2,Pmu,Wc,Wp,Yh,...
    Yalpha,Pa,Psi1,Wrho,Wd];

%% Setting up
% calculate the plastic work for isotropic compression path to reah the
% initial state
ConsolidationPressure = iConsolidationPressure+Tens;
% Initial stress level
isigma = isigma+Tens; sigma0 = isigma*[1,1,1,0,0,0]';
% flag works as regulator for hardening and softening regime
flag = 0; A = 0; B = 0;

[I1,I2,I3,~] = Invariant(sigma0);
iPlasticWork  = Wc*(I1/Pa).^Wp; PlasticWork0 = iPlasticWork; 
Q1 = zeros(1,consolidation_steps); P1 = zeros(1,consolidation_steps);
SIGMA = zeros(6,shearing_steps); Q2 = zeros(1,shearing_steps); P2 = zeros(1,shearing_steps);
CSTRAIN = zeros(6,consolidation_steps); SSTRAIN = zeros(6,shearing_steps);

%% isotropic compression
for i=1:consolidation_steps    
    [sigma1,PlasticWork1,DSTRAIN,flag,A,B] = SubsteppingStrainIncrementAlgorithm(dstrain1,PlasticWork0,m,b,consolidation_type,sigma0,flag,A,B,ModelParameter);
    sigma0 = sigma1; PlasticWork0 = PlasticWork1;
    [I1,~,~,J2] = Invariant(sigma0); q = sqrt(3*J2); p = I1;
    SIGMA(:,i) = sigma0; Q1(:,i) = q; P1(:,i) = p;
    if I1 >= 3*ConsolidationPressure
        break
    end
    CSTRAIN(:,i+1) = CSTRAIN(:,i)+DSTRAIN;
end

%% shearing process
for i=1:shearing_steps
    [sigma1,PlasticWork1,DSTRAIN,flag,A,B] = SubsteppingStrainIncrementAlgorithm(dstrain2,PlasticWork0,m,b,shearing_type,sigma0,flag,A,B,ModelParameter);
    sigma0 = sigma1; PlasticWork0 = PlasticWork1;
    [I1,~,~,J2] = Invariant(sigma0);
    % Octahedral shear stress
    q = sqrt(2*J2/3);
    p = I1; P2(:,i) = p;
    SIGMA(:,i) = sigma0; Q2(:,i+1) = q;
    SSTRAIN(:,i+1) = SSTRAIN(:,i) + DSTRAIN;
    CheckSigma = q/max(Q2);
    if q >= 10 && CheckSigma <= 0.6
        SSTRAIN(:,i+2:end) = [];
        Q2(:,i+2:end) = [];
        break
    end
end

%% Post processing
I1 = SSTRAIN(1,:)+SSTRAIN(2,:)+SSTRAIN(3,:);
J2 = (1/6)*((SSTRAIN(1,:)-SSTRAIN(2,:)).^2+(SSTRAIN(2,:)-SSTRAIN(3,:)).^2+(SSTRAIN(3,:)-SSTRAIN(1,:)).^2)+SSTRAIN(4,:).^2+SSTRAIN(5,:).^2+SSTRAIN(6,:).^2;

% Octahedral shear strain
GammaStrain = sqrt(8/3*J2); VolStrain = I1;
if shearing_type == 1
    filename = strcat('b = ',num2str(b),'.mat');
    save(filename)
else
    filename = strcat('sigma_c = ', num2str(iConsolidationPressure'),'.mat');
    save(filename)
end

figure('Color','w')
h1 = plot(GammaStrain,Q2,'k','LineWidth',4);
gca.LineWidth = 2.0;
xlabel('\gamma_{oct}')
ylabel('\tau_{oct}')
if shearing_type == 1
    title('Shear stress vs shear strain')
    legend(strcat('b = ',num2str(b)),'Location','SouthEast')
else
    title('Shear stress vs shear strain')
    legend(strcat('\sigma_c = ',num2str(iConsolidationPressure)),'Location','SouthEast')
end
grid on
saveas(h1,'Shear stress vs shear strain')

figure('Color','w')
h2 = plot(GammaStrain,VolStrain,'k','LineWidth',4);
gca.LineWidth = 2.0;
xlabel('\gamma_{oct}')
ylabel('\epsilon_{v}')
if shearing_type == 1
    title('Volumetric strain vs shear strain')
    legend(strcat('b = ',num2str(b)),'Location','SouthEast')
else
    title('Volumetric strain vs shear strain')
    legend(strcat('\sigma_c = ',num2str(iConsolidationPressure)),'Location','SouthEast')
end
grid on
saveas(h2,'Volumetric strain vs shear strain')
