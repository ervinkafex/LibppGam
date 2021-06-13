c written by Andrew M Taylor, Oct 2014
c ===========================================================
c     _     _ _                 _____
c    | |   (_) |               |  __ \
c    | |    _| |__  _ __  _ __ | |  \/ __ _ _ __ ___
c    | |   | | '_ \| '_ \| '_ \| | __ / _` | '_ ` _ \
c    | |___| | |_) | |_) | |_) | |_\ \ (_| | | | | | |
c    \_____/_|_.__/| .__/| .__/ \____/\__,_|_| |_| |_|
c                  | |   | |
c                  |_|   |_|
c ===========================================================
c If used please reference Phys.Rev. D90 (2014) 12, 123014 (astro-ph/1406.7369), "Parametrization of gamma-ray production cross-sections for pp interactions in a broad proton energy range from the kinematic threshold to PeV energies" by Ervin Kafexhiu, Felix Aharonian, Andrew M. Taylor, Gabriela S. Vila
c ===========================================================

      Program LibppGam

      implicit none
      Integer Models ! no. of hadronic models considered
      parameter (Models=4)
      Integer n,nmax,M,tracking,model
c nmax- photon energies
      parameter (nmax=1000) ! no. of photon energy bins
      Real*8 Eg(0:nmax),EgdNdE(1:Models,0:nmax),dEg
      Real*8 Eg2(0:nmax),Eg2dNdE(0:nmax)
      Real*8 a1(1:Models),a2(1:Models),a3(1:Models)
      Real*8 a4(1:Models),a5(1:Models)
      Real*8 b0,b1(1:Models),b2(1:Models),b3(1:Models)
      Real*8 factor,pi,K,s_inel
      Real*8 sigma1pi,sigma2pi,sigma0,sigmapi0,npi0
      Real*8 mpi,mp,s,eta,fBW,Mres,Gres,gamma
      Real*8 Tp,Tptilda,Tpth,Amax,M1,G1,Eg_min,Eg_max
      Real*8 Egmin,Egmax,Xgtilda,Xgtildam,Xg
      Real*8 Q,Q2,alpha,beta,C,C1,C2,C3,C4,B,mu
      Real*8 ds_dEg,Epim_CM,ppim_CM,Epim_LAB
      Character*90 outfile

      tracking=0

ccccccc input variable- proton KE, Tp [GeV] (to collide with proton at rest)
      Tp=1.0
cccccccccccccccccc

c min + max energy of output photon spectrum
      Eg_min=-5.0
      Eg_max=20.0

c zero photon spectra array for different hadronic models
      do M=1,Models
      do n=0,nmax
         EgdNdE(m,n)=0.0
      enddo
      enddo
   
      pi=4.0*atan(1.d0)

c pi-zero mass in GeV
      mpi=0.134976d0

c proton mass in GeV
      mp=0.938272d0

c resonance mass in GeV
      Mres=1.1883d0
      
c resonance width in GeV
      Gres=0.2264d0

      M1=0.25d0
      G1=0.15d0

      Tpth=2.0*mpi+mpi**2/(2.0*mp)

      gamma=sqrt(Mres**2*(Mres**2+Gres**2))

      K=sqrt(8.0)*Mres*Gres*gamma/(pi*sqrt(Mres**2+gamma))

c in mb
      sigma0=7.66e-3

      do n=0,nmax
         Eg(n)=10**(Eg_min+Eg_max*real(n)/real(nmax))
      enddo

      call get_parameters(Models,a1,a2,a3,a4,a5)

      write(outfile,'(A)')"EgdNdE.out"

c M=1 (Geant4), M=2 (Pythia 8), M=3 (SIBYLL), M=4 (QGSJET)
      do M=1,Models

c PHOTON ENERGIES
      do n=0,nmax

      dEg=Eg(n)*(10**(Eg_max/real(nmax))-1.0)

      call get_Eg_minmax(mpi,mp,Tp,s,Epim_LAB,Egmin,Egmax)

      if (tracking.eq.1) then
      print*,'Egmax:',Egmax,Tp-Tpth
      endif

      call get_parameters2(Tp,Models,b0,b1,b2,b3)

      call get_sigmapi0(Models,M,tracking,a1,a2,a3,a4,a5,
     -s,mpi,mp,K,Mres,Gres,sigma0,Tp,Tpth,s_inel,sigmapi0)

      Tptilda=Tp/mp

      call get_Amax(Models,M,mp,Epim_LAB,Tp,Tptilda,b0,b1,b2,b3,
     -sigmapi0,Amax)

      call get_parameters3(tracking,M,mp,mpi,Tp,Tptilda,Egmax,Q,Q2,
     -Xgtildam,C,C1,C2,C3,C4,mu,alpha,beta)

      Xgtilda=Eg(n)+(mpi**2.0)/(4.0*Eg(n))
      Xg=(Xgtilda-mpi)/(Xgtildam-mpi)

      call get_fX(M,Tp,ds_dEg,Amax,Xg,beta,alpha,C,C1,C2,C3,C4)      

      if (Tp.ge.Tpth.and.Eg(n).ge.Egmin.and.Eg(n).lt.Egmax) then
         if (isnan(ds_dEg)) print*,'Eg:',Eg(n),' Tp:',Tp
      endif

      if (Eg(n)*ds_dEg.le.1.0e-40) then
         EgdNdE(m,n)=EgdNdE(m,n)
      elseif (Tp.ge.Tpth.and.Eg(n).ge.Egmin.and.Eg(n).lt.Egmax) then
         EgdNdE(m,n)=EgdNdE(m,n)+Eg(n)*ds_dEg
      endif

c enddo to n do-loop (PHOTON ENERGY- Eg)
      enddo

c enddo to M do-loop (Models)
      enddo

      open(unit=10,file=outfile,status='unknown')

      write(10,'(A)')"# Eg [GeV],   EgdNdE (Geant4),   "
     -//"EgdNdE (pythia8),  EgdNdE (SIBYLL),   EgdNdE (QGSJET)"
      do n=0,nmax
 20      format(15(E10.4,9X))
         write(10,20)Eg(n),(EgdNdE(m,n),m=1,models)
      enddo
      close(10)
      
 99   end

cccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_parameters(Models,a1,a2,a3,a4,a5)

      implicit none
      Integer Models
      Real*8 a1(1:Models),a2(1:Models),a3(1:Models)
      Real*8 a4(1:Models),a5(1:Models)

      a1(1)=0.728
      a2(1)=0.596
      a3(1)=0.491
      a4(1)=0.2503
      a5(1)=0.117

      a1(2)=0.652
      a2(2)=0.0016
      a3(2)=0.488
      a4(2)=0.1928
      a5(2)=0.483

      a1(3)=5.436
      a2(3)=0.254
      a3(3)=0.072
      a4(3)=0.075
      a5(3)=0.166

      a1(4)=0.908
      a2(4)=9.0e-4
      a3(4)=6.089
      a4(4)=0.176
      a5(4)=0.448

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_fX(M,Tp,ds_dEg,Amax,Xg,beta,alpha,C,C1,C2,C3,C4)

      implicit none
      Integer M
      Real*8 Tp,ds_dEg,Amax,Xg,beta,alpha,C,C1,C2,C3,C4

c low energies
      if (Tp.lt.1.0) then
         ds_dEg=Amax*(1.0-Xg)**beta/((1.0+Xg/C)**(alpha))

c Geant4 (default)
      elseif (Tp.ge.1.0) then

         if (Tp.le.20.0) then
            ds_dEg=Amax*(1.0-Xg)**beta/((1.0+Xg/C1)**(alpha))
         elseif (Tp.gt.20.0) then
            ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C1)**(alpha))
         endif
         
      endif

c pythia8
      if (M.eq.2.and.Tp.gt.50.0) then
         ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C2)**(alpha))
c SIBYLL
      elseif (M.eq.3.and.Tp.gt.100.0) then
         ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C3)**(alpha))
c QGSJET
      elseif (M.eq.4.and.Tp.gt.100.0) then
         ds_dEg=Amax*(1.0-sqrt(Xg))**beta/((1.0+Xg/C4)**(alpha))
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_Eg_minmax(mpi,mp,Tp,s,Epim_LAB,Egmin,Egmax)

      implicit none
      Real*8 Tp,mpi,mp,s,Epim_CM,ppim_CM,gamma_CM,beta_CM
      Real*8 Epim_LAB,gamma_LAB,beta_LAB,Egmin,Egmax

      s=2.0*mp*(Tp+2.0*mp)

      Epim_CM=(s-4.0*mp**2+mpi**2)/(2.0*sqrt(s))
      ppim_CM=sqrt(Epim_CM**2-mpi**2)

      gamma_CM=(Tp+2.0*mp)/(sqrt(s))
      beta_CM=sqrt(1.0-gamma_CM**(-2))

      Epim_LAB=gamma_CM*(Epim_CM+beta_CM*ppim_CM)
      
      gamma_LAB=Epim_LAB/mpi
      beta_LAB=sqrt(1.0-gamma_LAB**(-2))

c      Egmax=Tp-Tpth ! approximate expression
      Egmin=(mpi/2.0)*gamma_LAB*(1.0-beta_LAB)
      Egmax=(mpi/2.0)*gamma_LAB*(1.0+beta_LAB)

      end

cccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_sigmapi0(Models,M,tracking,a1,a2,a3,a4,a5,
     -s,mpi,mp,K,Mres,Gres,sigma0,Tp,Tpth,s_inel,sigmapi0)
     
      implicit none
      Integer Models,M,tracking
      Real*8 eta,fBW,s,mpi,mp,K,Mres,Gres,Tp,Tpth
      Real*8 sigma0,sigma1pi,sigma2pi,s_inel,sigmapi0
      Real*8 Q,Q2,npi0
      Real*8 a1(1:Models),a2(1:Models),a3(1:Models)
      Real*8 a4(1:Models),a5(1:Models)

      eta=sqrt((s-mpi**2-4.0*mp**2)**2-16.0*mpi**2*mp**2)/
     -(2.0*mpi*sqrt(s))

c      fBW=K/(((sqrt(s)-mp)**2-Mres**2)**2+(Mres**2)*(Gres**2))
      fBW=mp*K/(((sqrt(s)-mp)**2-Mres**2)**2+(Mres**2)*(Gres**2))

      if (Tp.le.2.0) then
      sigma1pi=sigma0*eta**(1.95)*(1.0+eta+eta**5)*(fBW**1.86)
      endif

      s_inel=(30.7-0.96*log(Tp/Tpth)+0.18*(log(Tp/Tpth)**2))*
     -((1.0-(Tpth/Tp)**1.9)**3)

      if (Tp.lt.0.56) then
         sigma2pi=0.0
      elseif (Tp.ge.0.56.and.Tp.le.2.0) then
         sigma2pi=5.7/(1.0+exp(-9.3*(Tp-1.4)))
      endif

      if (tracking.eq.1) then
      print*,'s_inel:',s_inel
      endif

      if (Tp.lt.2.0) then
      sigmapi0=sigma1pi+sigma2pi
      elseif (Tp.ge.2.0.and.Tp.lt.5.0) then
      Q=(Tp-Tpth)/mp
      npi0=-6.0*10**(-3.0)+0.237*Q-0.023*Q**2
      sigmapi0=npi0*s_inel
      elseif (Tp.ge.5.0.and.Tp.le.100.0) then
      Q2=(Tp-3.0)/mp
      npi0=a1(1)*(Q2**a4(1))*(1.0+exp(-a2(1)*(Q2**a5(1))))*
     -(1.0-exp(-a3(1)*(Q2**0.25)))
      sigmapi0=npi0*s_inel
      elseif (Tp.gt.100.0) then
      Q2=(Tp-3.0)/mp
      npi0=a1(M)*(Q2**a4(M))*(1.0+exp(-a2(M)*(Q2**a5(M))))*
     -(1.0-exp(-a3(M)*(Q2**0.25)))
      sigmapi0=npi0*s_inel
      endif

      if (Tp.gt.50.0.and.M.eq.2) then
      Q2=(Tp-3.0)/mp
      npi0=a1(M)*(Q2**a4(M))*(1.0+exp(-a2(M)*(Q2**a5(M))))*
     -(1.0-exp(-a3(M)*(Q2**0.25)))
      sigmapi0=npi0*s_inel
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_Amax(Models,M,mp,Epim_LAB,Tp,Tptilda,b0,b1,b2,b3,
     -sigmapi0,Amax)

      implicit none
      Integer Models,M
      Real*8 mp,Epim_LAB,Tp,Tptilda,b0
      Real*8 b1(1:Models),b2(1:Models),b3(1:Models)
      Real*8 sigmapi0,Amax

c low energy maximum
      if (Tp.lt.1.0) then

      Amax=b0*sigmapi0/Epim_LAB

c Geant4 (default)
      elseif (Tp.ge.1.0.and.Tp.le.100.0) then

      Amax=b1(1)*Tptilda**(-b2(1))*exp(b3(1)*(log(Tptilda))**2)*
     -sigmapi0/mp      

      elseif (Tp.gt.100.0) then

      Amax=b1(M)*Tptilda**(-b2(M))*exp(b3(M)*(log(Tptilda))**2)*
     -sigmapi0/mp

      endif

c Pythia is exception since it can take over at 50 GeV
      if (Tp.gt.50.0.and.M.eq.2) then
      Amax=b1(M)*Tptilda**(-b2(M))*exp(b3(M)*(log(Tptilda))**2)*
     -sigmapi0/mp
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_parameters2(Tp,Models,b0,b1,b2,b3)

      implicit none
      Integer Models
      Real*8 Tp,b0,b1(1:Models),b2(1:Models),b3(1:Models)

      b0=5.9

      if (Tp.ge.1.0.and.Tp.lt.5.0) then
         b1(1)=9.53
         b2(1)=0.52
         b3(1)=0.054
      elseif (Tp.ge.5.0) then
         b1(1)=9.13
         b2(1)=0.35
         b3(1)=9.7e-3
      endif

      b1(2)=9.06
      b2(2)=0.3795
      b3(2)=0.01105
      
      b1(3)=10.77
      b2(3)=0.412
      b3(3)=0.01264

      b1(4)=13.16
      b2(4)=0.4419
      b3(4)=0.01439

      end

cccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine get_parameters3(tracking,M,mp,mpi,Tp,Tptilda,Egmax,
     -Q,Q2,Xgtildam,C,C1,C2,C3,C4,mu,alpha,beta)

      implicit none
      Integer tracking,M
      Real*8 mp,mpi,Tp,Tptilda,Egmax,Q,Q2,Xgtildam
      Real*8 C,C1,C2,C3,C4,mu,alpha,beta

      Q=(Tp-0.3)/mp
      Q2=(Tp-1.0)/mp
            
      Xgtildam=Egmax+(mpi**2)/(4.0*Egmax)

      if (tracking.eq.1) then
      print*,'Q:',Q
      endif

      C=5.0*mpi/Xgtildam
      C1=3.0*mpi/Xgtildam
      C2=3.5*mpi/Xgtildam
      C3=3.55*mpi/Xgtildam
      C4=3.55*mpi/Xgtildam

      mu=1.2*Q2**1.2*exp(-1.2*Q2)

      if (Tp.lt.1.0) then
         alpha=0.0
         beta=3.29-0.2*Tptilda**(-1.5)
c Geant4 (default)
      elseif (Tp.ge.1.0) then

      if (Tp.le.4.0) then
         alpha=mu+1.45
         beta=mu+2.45
      elseif (Tp.gt.4.0.and.Tp.le.20.0) then
         alpha=mu+1.5
         beta=1.5*mu+4.95
      elseif (Tp.gt.20.0.and.Tp.le.100.0) then
         alpha=1.0
c         beta=4.4
         beta=4.2
      elseif (Tp.gt.100.0) then
         alpha=1.0
         beta=4.9
      endif

      endif

c pythia8
      if (M.eq.2.and.Tp.gt.50.0) then
         alpha=1.0
         beta=4.00001
c SIBYLL
      elseif (M.eq.3.and.Tp.gt.100.0) then
         alpha=1.0
         beta=3.6
c QGSJET
      elseif (M.eq.4.and.Tp.gt.100.0) then
         alpha=1.0
         beta=4.5
      endif

      if (tracking.eq.1) then
      print*,'C:',C,'Xgtildam:',Xgtildam
      print*,'alpha:',alpha,' beta:',beta
      endif

      end
