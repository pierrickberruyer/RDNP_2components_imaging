% Berruyer, Pierrick* ; Cibaka-Ndaya, Cynthia ; Pinon, Arthur ; Sanchez, Clément ; Drisko, Glenna L. ; Emsley, Lyndon*
% * = corresponding authors
% Available under the license CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International). See licensing details at: https://creativecommons.org/licenses/by-sa/4.0/. 
% In case of use, please cite: DOI: 10.5281/zenodo.7756308 and DOI: 10.1021/jacs.3c01279
% Disclaimer: please note the code is provided free of charge and come with no warranty. 

% Based on the matlab codes published in 
% Rossini, A. J.; Zagdoun, A.; Hegner, F.; Schwarzwälder, M.; Gajan, D.; Copéret, C.; Lesage, A.; Emsley, L. Dynamic Nuclear Polarization NMR Spectroscopy of Microcrystalline Solids. J. Am. Chem. Soc. 2012, 134 (40), 16899–16908. https://doi.org/10.1021/ja308135r.
% Pinon, A. C.; Schlagnitweit, J.; Berruyer, P.; Rossini, A. J.; Lelli, M.; Socie, E.; Tang, M.; Pham, T.; Lesage, A.; Schantz, S.; et al. Measuring Nano- to Microstructures from Relayed Dynamic Nuclear Polarization NMR. J. Phys. Chem. C 2017, 121 (29), 15993–16005. https://doi.org/10.1021/acs.jpcc.7b04438.
% Berruyer, P.; Björgvinsdóttir, S.; Bertarello, A.; Stevanato, G.; Rao, Y.; Karthikeyan, G.; Casano, G.; Ouari, O.; Lelli, M.; Reiter, C.; et al. Dynamic Nuclear Polarization Enhancement of 200 at 21.15 T Enabled by 65 KHz Magic Angle Spinning. J. Phys. Chem. Lett. 2020, 11 (19), 8386–8391. https://doi.org/10.1021/acs.jpclett.0c02493.

clear
clc

%% Loading DNP enhancement build-up data
%% Experimental data for 2 components PS/SiO2-CTAB (PS=comp1 and CTAB = comp2) core-shell particles can be obtained from 10.5281/zenodo.7756308 and DOI: 10.1021/jacs.3c01279
%% enhancements  matrix(3:length(time)) organised as (juice ; _comp1 ; _comp2)x(time)
if exist('SampleA_2_extracted_nmr_data.mat','file')
    load('SampleA_2_extracted_nmr_data.mat','time','enhancements')
else
    time = linspace(0,20,12)'; enhancements = rand(3,12);
    warning('Using synthetic random data')
end

Enhancement_comp2 = enhancements(3,1:12);
Enhancement_comp1 = enhancements(2,1:12);
EnhancementJuice = enhancements(1,1:12) ;
delays = time(1:12);
clear enhancements time

%% Create parallel computing pool 
if isempty(gcp('nocreate'))
    thepool = parpool('local');
end

%% Setting max time of calculation before forcing to stop if not converge
maxtimecluster = 24*3600 - 5*60; % max time in seconds (minus 5 min for safety - matlab and parpool starting time)


%% Define number of total pixel used in the simulation (particle + juice layer)
totalpixels = 200;

%% Define number of pixel corresponding to the particle (impose the radius of the full particle)
totalpixelsmaterial = 80 ; % particle radius size = 80 nm for sample A

%% Define all constant RDNP parameters
m = 2; % spherical symmetry
totallength = 0.2; % total length of the simulation in um (particle + juice)
rstep = 400; % number of points in the space domain
T1juice = 2.5; % Juice TB in s
T1_comp2 = 52; % _comp2 T1 in s
T1_comp1 = 1.3; % _comp1 T1 in s
epsilon0 = 140; % Juice DNP enhancement
theta = 1; % Percentage of QUENCHED signal in the source (0 if no quenched signal)
dep = 0.4; % Percentage of depolarization in the source
concjuice = 44; % 1H concentration in the source
conc_comp2 = 57/2; % 1H concentration of _comp2 in the layer (pure _comp2 = 57 M)
conc_comp1 = 77; %1H concentration in _comp1
MAS = 8; % MAS rate in kHz
D0 = 4.25e-3; % Ref diffusion coefficient of _comp1 in a static sample
D_comp1 = D0*(conc_comp1/70)^(1/3)/(1+0.6*MAS); % The diffusion rate is scaled with the proton concentration and MAS
D_comp2 = D0*(conc_comp2/70)^(1/3)/(1+0.6*MAS); % The diffusion rate is scaled with the proton concentration and MAS
Djuice = D0*(concjuice/70)^(1/3)/(1+0.6*MAS); % The diffusion rate is scaled with the proton concentration and MAS
Cp = 7.98e-2; % Molar heat capacity of 1H spins (should not be changed unless the field is different from 9.4 T)

%% Generate a random composition for the particle as the starting guess composition
simplematerial_start = [rand(1,totalpixelsmaterial)]; % (1 x totalpixelsmaterial) matrix

%% Define upper and lower boundaries of the composition function (ie. between 0 and 1 by definition of the C_k(r)
discretcomposition_lb = zeros(size(simplematerial_start)); 
discretcomposition_ub = ones(size(simplematerial_start));

%% Calculate initial residue
sse1 = sse_imaging(simplematerial_start,Enhancement_comp1,Enhancement_comp2,EnhancementJuice,delays,m,totallength,rstep,T1juice,T1_comp2,T1_comp1,epsilon0,theta,dep,concjuice,conc_comp2,conc_comp1,D_comp1,D_comp2,Djuice,Cp,totalpixels,totalpixelsmaterial)

%% Imaging the particle internal composition using Pattern Search
sse_imagingfun = @(discretcomposition)sse_imaging(discretcomposition,Enhancement_comp1,Enhancement_comp2,EnhancementJuice,delays,m,totallength,rstep,T1juice,T1_comp2,T1_comp1,epsilon0,theta,dep,concjuice,conc_comp2,conc_comp1,D_comp1,D_comp2,Djuice,Cp,totalpixels,totalpixelsmaterial);
options = optimoptions('patternsearch','Display','iter','UseParallel',true,'InitialMeshSize',1,'MaxMeshSize',1,'MaxTime',maxtimecluster,'MaxFunctionEvaluations',100000000000000,'MaxIterations',60000);
fit_result = patternsearch(sse_imagingfun,simplematerial_start,[],[],[],[],discretcomposition_lb,discretcomposition_ub,[],options);

%% Saving the whole Matlab final space
save final_matlab_workspace_composition_fitted
fprintf('Final results saved!\n');


%% Definition of functions

function [sse] = sse_imaging(discretcomposition,Enhancement_comp1,Enhancement_comp2,EnhancementJuice,delays,m,totallength,rstep,T1juice,T1_comp2,T1_comp1,epsilon0,theta,dep,concjuice,conc_comp2,conc_comp1,D_comp1,D_comp2,Djuice,Cp,totalpixels,totalpixelsmaterial)
% Imaging residue function
% Compute the residue between the simulated enhancements as function of time for a
% given composition and the experimental enhancements.

totalfixedjuicepixel = totalpixels - totalpixelsmaterial;

% Creating the matrix of composition from the simulation input where the
% composition is given in matrix (1 x 3*totalpixelsmaterial)  and transform
% to (3 x totalpixelsmaterial) ordered as   _comp1
%                                           _comp2
%                                           Juice
discretcomposition3d = abs([discretcomposition ; ones(1,totalpixelsmaterial)-discretcomposition ; zeros(1,totalpixelsmaterial)]);

% Adding to the composition a layer of pure juice
constantjuicelayer = [zeros(1,totalfixedjuicepixel) ; zeros(1,totalfixedjuicepixel) ; ones(1,totalfixedjuicepixel)];
discretcomposition3d = [discretcomposition3d constantjuicelayer];


[epsilon_comp1,epsilon_comp2,epsilonJuice] = RDNPfit(m,rstep,delays,T1juice,T1_comp2,T1_comp1,epsilon0,D_comp1,D_comp2,Djuice,theta,dep,concjuice,conc_comp2,conc_comp1,Cp,totallength,discretcomposition3d);
 

   ErrorVector_comp1 =  epsilon_comp1(2:11)  - Enhancement_comp1(2:11);
   ErrorVector_comp2 =  epsilon_comp2(2:11)  - Enhancement_comp2(2:11);
   ErrorVectorJuice =  epsilonJuice(2:11)  - EnhancementJuice(2:11);
    
    s_comp1 = sum(ErrorVector_comp1.^ 2);
    s_comp2 = sum(ErrorVector_comp2.^ 2);
    sJuice = sum(ErrorVectorJuice.^ 2);

    sse = sqrt(s_comp1 + s_comp2 + sJuice) ; 
    
end

function [epsilon_comp1,epsilon_comp2,epsilonJuice]=RDNPfit(m,rstep,delays,T1juice,T1_comp2,T1_comp1,epsilon0,D_comp1,D_comp2,Djuice,theta,dep,concjuice,conc_comp2,conc_comp1,Cp,totallength,discretcomposition)

 r = linspace(0,totallength,rstep); 
 t = delays;

 
% Extract the first solution of the pdepe equation, as p and pref.
 sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,r,t); %solution of the diffusion equation with microwaves
   p = sol(:,:,1); %p is the polarization with microwaves
sol2 = pdepe(m,@pdex1pdenoDNP,@pdex1icnoDNP,@pdex1bcnoDNP,r,t);%solution of the diffusion equation without microwaves
pref = sol2(:,:,1); %pref is the polarization without microwaves

% Process the polarization and calculate the enhancement
  
    %Initialization of matrixes
    
    p_cut = zeros(length(t),length(r));
    pref_cut = zeros(length(t),length(r));
    Presign = zeros(length(t),length(r));
    Presignref = zeros(length(t),length(r));
    Presign_comp1 = zeros(length(t),length(r));
    Presignref_comp1 = zeros(length(t),length(r));
    Presign_comp2 = zeros(length(t),length(r));
    Presignref_comp2 = zeros(length(t),length(r));
    PresignJuice = zeros(length(t),length(r));
    PresignrefJuice = zeros(length(t),length(r));    
    SignalPol = zeros(1,length(t));
    SignalRef = zeros(1,length(t));
    SignalPol_comp1 = zeros(1,length(t));
    SignalRef_comp1 = zeros(1,length(t));
    SignalPol_comp2 = zeros(1,length(t));
    SignalRef_comp2 = zeros(1,length(t));
    SignalPolJuice = zeros(1,length(t));
    SignalRefJuice = zeros(1,length(t));    
    NormSignalPol = zeros(1,length(t));
    NormSignalRef = zeros(1,length(t));
    NormSignalPol_comp1 = zeros(1,length(t));
    NormSignalRef_comp1 = zeros(1,length(t)); 
    NormSignalPol_comp2 = zeros(1,length(t));
    NormSignalRef_comp2 = zeros(1,length(t));
    NormSignalPolJuice = zeros(1,length(t));    
    NormSignalRefJuice = zeros(1,length(t));
    epsilon = zeros(1,length(t));
    epsilon_comp1 = zeros(1,length(t));   
    epsilon_comp2 = zeros(1,length(t));
    epsilonJuice = zeros(1,length(t));
    Normepsilon = zeros(1,length(t));


    for l=1:length(r)
        p_cut(:,l) = p(:,l)*theta0(r(l)); %Visible polarization ON
        pref_cut(:,l) = pref(:,l)*theta0(r(l));  %Visible polarization OFF
    end



% Calcul of Presignal
    for j=1:length(t)
        for k=1:length(r) % EVERYTHING     
                Presign(j,k) = fc(r(k))*p_cut(j,k)*(r(k))^2; %3D
                Presignref(j,k) = fc(r(k))*pref_cut(j,k)*(r(k))^2; 
        end  
    end
    
    for j=1:length(t)
        for k=1:length(r) % _comp1 
                Presign_comp1(j,k) = fc_comp1(r(k))*p_cut(j,k)*(r(k))^2; %3D
                Presignref_comp1(j,k) = fc_comp1(r(k))*pref_cut(j,k)*(r(k))^2; 
        end  
    end   
    for j=1:length(t)
        for k=1:length(r) % _comp2   
                Presign_comp2(j,k) = fc_comp2(r(k))*p_cut(j,k)*(r(k))^2; %3D
                Presignref_comp2(j,k) = fc_comp2(r(k))*pref_cut(j,k)*(r(k))^2; 
        end  
    end
 
    for j=1:length(t)
        for k=1:length(r) % Juice   
                PresignJuice(j,k) = fcJuice(r(k))*p_cut(j,k)*(r(k))^2; %3D
                PresignrefJuice(j,k) = fcJuice(r(k))*pref_cut(j,k)*(r(k))^2; 
        end  
    end
    
   
% Calcul of the Signals
    for w=1:length(t)
    SignalPol(w) = trapz(Presign(w,:));
    SignalRef(w) = trapz(Presignref(w,:));
    SignalPol_comp1(w) = trapz(Presign_comp1(w,:));
    SignalRef_comp1(w) = trapz(Presignref_comp1(w,:));
    SignalPol_comp2(w) = trapz(Presign_comp2(w,:));
    SignalRef_comp2(w) = trapz(Presignref_comp2(w,:));
    SignalPolJuice(w) = trapz(PresignJuice(w,:));
    SignalRefJuice(w) = trapz(PresignrefJuice(w,:));
    end 
    
% The normalized values of Signal Pol and Ref 
    for c=1:length(t)
    NormSignalPol(c) = SignalPol(c)./SignalPol(end); 
    NormSignalRef(c) = SignalRef(c)./SignalRef(end);
    NormSignalPol_comp1(c) = SignalPol_comp1(c)./SignalPol_comp1(end); 
    NormSignalRef_comp1(c) = SignalRef_comp1(c)./SignalRef_comp1(end);
    NormSignalPol_comp2(c) = SignalPol_comp2(c)./SignalPol_comp2(end); 
    NormSignalRef_comp2(c) = SignalRef_comp2(c)./SignalRef_comp2(end);
    NormSignalPolJuice(c) = SignalPolJuice(c)./SignalPolJuice(end); 
    NormSignalRefJuice(c) = SignalRefJuice(c)./SignalRefJuice(end);
    end
   
% The ratio of Pol and Ref gives the enhancement
 
    for v=1:length(t)
    epsilon(v) = SignalPol(v)./SignalRef(v);
    Normepsilon(v) = epsilon(v)./epsilon(end);
    epsilon_comp1(v) = SignalPol_comp1(v)./SignalRef_comp1(v);
    epsilon_comp2(v) = SignalPol_comp2(v)./SignalRef_comp2(v);
    epsilonJuice(v) = SignalPolJuice(v)./SignalRefJuice(v);
    end
   

% -------------------------------------------------------------------------
% Function definitions
% -------------------------------------------------------------------------

% Source term WITH DNP.
function [c,f,s] = pdex1pde(x,~,u,DuDx)
    c = Cp*fc(x);
    f = fD(x)*Cp*fc(x)*DuDx;
    s = -(u-fe0(x))*(1/fT1(x))*Cp*fc(x);
end
 
% Source term WITHOUT DNP.
function [c,f,s] = pdex1pdenoDNP(x,~,u,DuDx)
    c = Cp*fc(x);
    f = fD(x)*Cp*fc(x)*DuDx;
    s = -(u-depo(x))*(1/fT1(x))*Cp*fc(x);
end
 
% Initial polarization is 0 all over the crystal WITH DNP
function u0 = pdex1ic(~)
    u0 = 0;
end
 
% Initial polarization is 0 all over the crystal WITHOUT DNP
function u0 = pdex1icnoDNP(~)
    u0 = 0;
end
    
% Boundary conditions WITH DNP
function [pl,ql,pr,qr] = pdex1bc(~,~,~,~,~)  
   pl = 0; 
   ql = 1;
   pr = 0;
   qr = 1;
end
    
% Boundary condition WITHOUT DNP
   function [pl,ql,pr,qr] = pdex1bcnoDNP(~,~,~,~,~)  
   pl = 0;
   ql = 1;
   pr = 0;
   qr = 1;
   end
   
% T1 function that will determine the relaxation term
function [T1] = fT1(y)
    [comp1,comp2,juice] = fcomposition(y,totallength,discretcomposition);
    T1 = T1_comp1*comp1 + T1_comp2*comp2 + T1juice*juice;
end


% Build-up epsilon function MW ON
 function [epsilonzero] = fe0(y)
        [~,~,juice] = fcomposition(y,totallength,discretcomposition);
        epsilonzero = 1+ ((epsilon0-1)*dep)*juice;
 end
        
% Percentage of seen signal
 function [theta3] = theta0(y)
        [~,~,juice] = fcomposition(y,totallength,discretcomposition);
        theta3 = 1 + theta*juice;    

 end    
% D function that depends on the concentration of protons
function [D] = fD(y)
         [comp1,comp2,juice] = fcomposition(y,totallength,discretcomposition);
         D = D_comp1*comp1 + D_comp2*comp2 + Djuice*juice ;
end
        
% Build-up function MW OFF (depolarization)
 function [depol] = depo(y)
        [~,~,juice] = fcomposition(y,totallength,discretcomposition);
        depol = 1 - (1-dep)*juice;    
 end

% Concentration spatial function
function [ctot] = fc(y)
    [comp1,comp2,juice] = fcomposition(y,totallength,discretcomposition);
    ctot = comp1*conc_comp1 + comp2*conc_comp2 + juice*concjuice;
        
end

function [c_comp1] = fc_comp1(y)
    [comp1,~,~] = fcomposition(y,totallength,discretcomposition);
    c_comp1 = comp1*conc_comp1;
        
end

function [c_comp2] = fc_comp2(y)
    [~,comp2,~] = fcomposition(y,totallength,discretcomposition);
    c_comp2 = comp2*conc_comp2;
        
end

function [cJuice] = fcJuice(y)
    [~,~,juice] = fcomposition(y,totallength,discretcomposition);
    cJuice = juice*concjuice;
        
end

end

function [comp1,comp2,juice] = fcomposition(y,totallength,discretcomposition)
%Computed a smooth function of composition based on a discret generated composition of comp1 in a material of lengh totallength.
% lengthmat is a scaller
% discretcomposition_comp1 is a 3 x nbpixels matrix where the first line give
% the comp1 
%, the second give the comp2 
%, and the 3rd line give the juice 
% y is the space variable

slopeparam = 3; 

nbpixels = size(discretcomposition,2);
sizepixels = totallength/nbpixels;
slope = slopeparam/(sizepixels);

plouf_comp1 = 0 ;
plouf_comp2 = 0 ;
ploufJuice = 0 ;

for i = 1:nbpixels
    
    if i == 1 
    centrepx = - sizepixels/2 ;    
    plouf_comp1 = plouf_comp1 + 0.5.*discretcomposition(1,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));
    plouf_comp2 = plouf_comp2 + 0.5.*discretcomposition(2,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));
    ploufJuice = ploufJuice + 0.5.*discretcomposition(3,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));    
    end
    
    centrepx = i*sizepixels - sizepixels/2 ;
    plouf_comp1 = plouf_comp1 + 0.5.*discretcomposition(1,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));
    plouf_comp2 = plouf_comp2 + 0.5.*discretcomposition(2,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));
    ploufJuice = ploufJuice + 0.5.*discretcomposition(3,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));

    if i == nbpixels 
     centrepx = i*sizepixels + sizepixels/2 ;   
    plouf_comp1 = plouf_comp1 + 0.5.*discretcomposition(1,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));
    plouf_comp2 = plouf_comp2 + 0.5.*discretcomposition(2,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));
    ploufJuice = ploufJuice + 0.5.*discretcomposition(3,i).*(tanh(slope*(y-(centrepx - sizepixels/2))) - tanh(slope*(y-(centrepx + sizepixels/2))));    
    end
    
    
end

comp1 = plouf_comp1;
comp2 = plouf_comp2 ;
juice = ploufJuice ;

end
