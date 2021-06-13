function F=multip_pi0_Pythia8(Tp)
% pi0 multiplicity from Geant4 and Pythia8.
% Tp is in GeV
% 
mpi = 0.134976; % GeV (pi0 mass)
mp  = 0.938272; % GeV (proton mass)
Tpth= 2*mpi + mpi^2/mp/2;
% ++++++++++++++++++++++++++++++++++
 F = zeros(size(Tp));
 
 ilg4=Tp<5;
 ihg4=(Tp>=5) & (Tp<=50);
 ih  =Tp>50;
% ++++++++++++++++++++++++++++++++++
% 1<=Tp<5 GeV Geant4
Qlg4 = (Tp(ilg4) - Tpth)/mp;

ALg4 = [-6e-3,0.237,-0.023];

F(ilg4) = ALg4(1) + ALg4(2)*Qlg4 + ALg4(3)*Qlg4.^2;

% ++++++++++++++++++++++++++++++++++
% 5<=Tp<=50 GeV Geant4
% A=[1.0917,1.2538,0.5312,0.2236];

AHg4=[0.728,0.596,0.491,0.2503,0.117];

Qhg4  = (Tp(ihg4)-3)/mp;

F(ihg4)= AHg4(1).*Qhg4.^AHg4(4).*...
(1 + exp(-AHg4(2)*Qhg4.^AHg4(5))).*(1-exp(-AHg4(3)*Qhg4.^0.25));
% ++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++
% Tp>50 GeV Pythia8

AH=[0.652,0.0016,0.488,0.1928,0.483];

Qh  = (Tp(ih)-3)/mp;

F(ih)= AH(1).*Qh.^AH(4).*...
(1 + exp(-AH(2)*Qh.^AH(5))).*(1-exp(-AH(3)*Qh.^0.25));

% ++++++++++++++++++++++++++++++++++

end
