function Amax=Amax_Geant4(Tp)
% gamma spectra Amax from Geant4.
% Tp is in GeV
% 
mpi = 0.134976; % GeV (pi0 mass)
mp  = 0.938272; % GeV (proton mass)
% ++++++++++++++++++++++++++++++++++
 Amax = zeros(size(Tp));
 
 iexp = Tp<1;
 ilg4 = (Tp>=1) & (Tp<5);
 ihg4 = (Tp>=5);
 
Tpmp  = Tp/mp;
lTpmp = log(Tpmp); 
 
% ++++++++++++++++++++++++++++++++++
% Tp<1 GeV From experiemental data fit
       ss = 2*mp*(Tp(iexp)+2*mp); %GeV^2
 EpimaxCM = (ss -4*mp^2 + mpi^2)./(2*sqrt(ss)); %GeV
 PpimaxCM = sqrt(EpimaxCM.^2 - mpi^2); %GeV
  gammaCM = (Tp(iexp) + 2*mp)./sqrt(ss);
   betaCM = sqrt(1 - 1./gammaCM.^2);
EpimaxLAB = gammaCM.*( EpimaxCM + betaCM.*PpimaxCM); %GeV

b0 = 5.9;

Amax(iexp) = b0 *mp./EpimaxLAB;
 
% ++++++++++++++++++++++++++++++++++
% 1<=Tp<5 GeV from Geant4 fit

ALg4 = [0.054,-0.52,9.53];
Amax(ilg4) = ALg4(3)*Tpmp(ilg4).^ALg4(2).*exp(ALg4(1)*lTpmp(ilg4).^2);

% ++++++++++++++++++++++++++++++++++
% 5<=Tp<=1e5 GeV from Geant4 fit
AHg4=[9.7e-3,-0.35,9.13];
Amax(ihg4) = AHg4(3)*Tpmp(ihg4).^AHg4(2).*exp(AHg4(1)*lTpmp(ihg4).^2);

% ++++++++++++++++++++++++++++++++++

Amax = Amax.*XSInlcPi0Geant4(Tp)/mp; %[GeV^-1]


end
