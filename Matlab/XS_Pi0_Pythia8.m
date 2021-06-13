function xs=XS_Pi0_Pythia8(Tp)
% This function calculates the pi0 production cross section.
% This cross section is constructed as follows:
% for energies Tp<=2GeV, includes fitting from 
% experimental data of the following channels:
% 1) p+p->p+p+pi0
% 2) p+p->p+n+(pi+)+pi0
% 3) p+p->p+p+2pi0
% 
% The fitting is valid from threshold to Tp=2 GeV.
% 
% From 2 GeV < Tp <= 50 GeV the cross section is taken from 
% Geant4 pi0 production multiplicity. For Tp > 50 GeV we take
% from Pythia 8 multiplicity.
% The cross section is in mb and proton kinetic energy in GeV.
% ==============================================================

mpi = 0.134976; % GeV (pi0 mass)
mp  = 0.938272; % GeV (proton mass)
Tpth= 2*mpi + mpi^2/mp/2;
% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo

lTp = size(Tp);

xs1pi     = zeros(lTp);
xs2pi     = zeros(lTp);
XS5toHigh = zeros(lTp);

% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo
% p+p -> p+p+ pi0
% One pion production from threshold to few GeV

M   = 1.1883; % resonance effective mass in GeV
Gam = 0.2264;% resonance effective width in GeV

% relativistic Breit-Wigner distribution
g1 = sqrt(M.^2.*(M.^2 + Gam.^2));
k  = sqrt(8)*M*Gam*g1/(pi*sqrt(M^2+g1));

iTpth2 = find(Tp>Tpth & Tp<=2);

ss  = 2*mp*(Tp(iTpth2) + 2*mp); % GeV^2
eta = sqrt((ss -mpi^2 -4*mp^2).^2 - 16*mpi^2*mp^2)./(2*mpi*sqrt(ss));
x   = sqrt(ss)-mp;

fBW = mp*k./((x.^2 - M^2).^2 + M^2*Gam^2);

sig0 = 7.66e-3;% mb

xs1pi(iTpth2)=sig0*eta.^1.95.*(1 + eta + eta.^5).*fBW.^1.86;% in mb

% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo
% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo
% Two pion event from threshold to few GeV
% % p+p -> p+n +(pi+) + pi0 & p+p -> D + (pi+) + pi0
% % p+p -> p+p + 2 pi0 & p+p -> p+p + pi0

i2pi=find(Tp>=0.56 & Tp<=2);

xs2pi(i2pi) = 5.7./(1+exp(-9.3*(Tp(i2pi)-1.4))); % in mb

% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo
% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo
% From Pythia8 multiplicity Tp > 2 GeV

i3pi = Tp>2;

XS5toHigh(i3pi) = Geant4Pythia8Pi0Multiplicity(Tp(i3pi)).*XSinel(Tp(i3pi));

% oooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOoooooooooooooooo

xs = xs1pi + xs2pi + XS5toHigh; % in mb

end
