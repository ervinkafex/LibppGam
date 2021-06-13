function dXSdEg=dXSdEg_SIBYLL(Tp,Egam)

dXSdEg = AmaxSIBYLL(Tp).* F_Shape_SIBYLL(Tp,Egam);

end

% =====================================================

function FF = F_Shape_SIBYLL(Tp,Egam)
mpi = 0.134976; % GeV (pi0 mass)
mp  = 0.938272; % GeV (proton mass) 
% --------------
% =====================================================
     Tpth = 2*mpi + mpi^2/mp/2;
       ss = 2*mp*(Tp+2*mp); %GeV^2
 EpimaxCM = (ss -4*mp^2 + mpi^2)./(2*sqrt(ss)); %GeV
 PpimaxCM = sqrt(EpimaxCM.^2 - mpi^2); %GeV
  gammaCM = (Tp + 2*mp)./sqrt(ss);
   betaCM = sqrt(1 - 1./gammaCM.^2);
EpimaxLAB = gammaCM.*( EpimaxCM + betaCM.*PpimaxCM); %GeV
gam_piLAB = EpimaxLAB/mpi;
bet_piLAB = sqrt(1-1./gam_piLAB.^2);
  Egammax = mpi/2* gam_piLAB.*(1+bet_piLAB);
%   Egammin = mpi/2* gam_piLAB.*(1-bet_piLAB);
% =====================================================

   A0 = zeros(size(Tp));
   FF = A0;
 beta = A0;
gamma = A0;
kappa1= A0;
XX1=A0;XX2=A0;XX3=A0;XX4=A0;XX5=A0;
B1 = ones(size(Tp));
CC2=B1;CC3=B1;CC4=B1;CC5=B1;
% =====================================================
Y  = Egam + mpi^2./Egam/4;
Y0 = Egammax + mpi^2./Egammax/4;
X  = (Y-mpi)./(Y0-mpi);
% =====================================================

theta = Tp/mp;
kappa = 3.29 - 0.2*theta.^(-1.5);
% --------------
q  = (Tp -1.0)/mp;
mu = 1.25*q.^(1.25).*exp(-1.25*q);
Cg4= 3*mpi./Y0;
C  = 3.55*mpi./Y0;
% =====================================================

i1 = (X<1) & (Tp>Tpth) & (Tp<1);
i2 = (X<1) & (Tp>=1) & (Tp<=4);
i3 = (X<1) & (Tp> 4) & (Tp<=20);
i4 = (X<1) & (Tp>20) & (Tp<=100);
i5 = (X<1) & (Tp>100);

% ------------------------------
XX1(i1)= X(i1);
kappa1(i1) = kappa(i1);
FF1 = (1-XX1).^kappa1;
% ------------------------------
 XX2(i2) = X(i2);
 CC2(i2) = Cg4(i2);
beta(i2) = mu(i2) + 2.45;
gamma(i2)= mu(i2) + 1.45;
FF2 = (1-XX2).^beta./(1+XX2./CC2).^gamma;
% ------------------------------
 XX3(i3) = X(i3);
 CC3(i3) = Cg4(i3);
beta(i3) = 1.5*mu(i3) + 4.95;
gamma(i3)= mu(i3) + 1.5;
FF3 = (1-XX3).^beta./(1+XX3./CC3).^gamma;
% ------------------------------
XX4(i4)= X(i4);
CC4(i4)= Cg4(i4);
FF4 = (1-sqrt(XX4)).^4.2./(1+XX4./CC4);   
% ------------------------------
XX5(i5)= X(i5);
CC5(i5)= C(i5);
FF5 = (1-sqrt(XX5)).^3.6./(1+XX5./CC5);   
% ------------------------------

FF(i1) = FF1(i1);
FF(i2) = FF2(i2);
FF(i3) = FF3(i3); 
FF(i4) = FF4(i4);
FF(i5) = FF5(i5);

end
