function F=multip_pi0_Geant4(Tp)
% pi0 multiplicity from Geant4.
% Tp is in GeV
% 
mpi = 0.134976; % GeV (pi0 mass)
mp  = 0.938272; % GeV (proton mass)
Tpth= 2*mpi + mpi^2/mp/2;
% ++++++++++++++++++++++++++++++++++
 F = zeros(size(Tp));
 
 il=Tp<5;
 ih=Tp>=5;
% ++++++++++++++++++++++++++++++++++
% 1<=Tp<5 GeV
Ql = (Tp(il) - Tpth)/mp;

AL = [-6e-3,0.237,-0.023];

F(il) = AL(1) + AL(2)*Ql + AL(3)*Ql.^2;

% ++++++++++++++++++++++++++++++++++
% 5<=Tp<=1e5 GeV
% A=[1.0917,1.2538,0.5312,0.2236];

AH=[0.728,0.596,0.491,0.2503,0.117];

Qh  = (Tp(ih)-3)/mp;

F(ih)= AH(1).*Qh.^AH(4).*...
(1 + exp(-AH(2)*Qh.^AH(5))).*(1-exp(-AH(3)*Qh.^0.25));
% ++++++++++++++++++++++++++++++++++


end
