function xs =XS_inel(Tp)
% 
% This function calculates the total inelastic cross section
% from our fit of the PDG data including TOTEM@LHC
% 
% Tp - is proton kinetic energy in the LAB frame in [GeV]
% xs - is the inelastic cross section in [mb]
% ===========================================================

mpi = 0.134976; % GeV (pi0 mass)
mp  = 0.938272; % GeV (proton mass)

Tpth = 2*mpi + mpi^2/(2*mp);

yy = Tp/Tpth;
x  = log(yy);

xs = (30.7 - 0.96*x + 0.18*x.^2).*((max(0.,1-1./yy.^1.9)).^3);

end
