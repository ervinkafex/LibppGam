/*
Writen by Ervin Kafexhiu, Oct 2014, and updated in May 2021
=====================================================================
     _     _ _                 _____
    | |   (_) |               |  __ \
    | |    _| |__  _ __  _ __ | |  \/ __ _ _ __ ___
    | |   | | '_ \| '_ \| '_ \| | __ / _` | '_ ` _ \
    | |___| | |_) | |_) | |_) | |_\ \ (_| | | | | | |
    \_____/_|_.__/| .__/| .__/ \____/\__,_|_| |_| |_|
                  | |   | |
                  |_|   |_|
=====================================================================
If used please reference to Phys.Rev. D90 (2014) 12, 123014 
(astro-ph/1406.7369), "Parametrization of gamma-ray 
production cross-sections for pp interactions in a broad 
proton energy range from the kinematic threshold to PeV 
energies" by 
Ervin Kafexhiu, Felix Aharonian, Andrew M. Taylor, Gabriela S. Vila
=====================================================================
This library contains the parametrization for different quantities 
for pp->pi0->2gamma. The basic units in each functions are GeV and mb.
*/

#ifndef LIBPPGAM_H
#define LIBPPGAM_H

//=======================================================
// Tp - is the proton kinetic energy in [GeV].
// Egamma - is the gamma-ray energy in [GeV].
//=======================================================

#define pi 3.141592653589793 //number pi.
#define m_p 0.938272 // GeV, proton mass (taken from PDG)
#define m_pi 0.134976 // GeV, pi0 mass (taken from PDG)

// the proton threshold energy in the LAB frame
#define Tp_th (2.0*m_pi + m_pi*m_pi/(2.0*m_p)) // in GeV


//=======================================================
// Calculates the maximum pi0 energy in the LAB frame.
double Epi0_max_LAB(double Tp);
// Calculates the maximum gamma-ray energy in the LAB frame.
double Egamma_max(double Tp);
//=======================================================
// pp total inelastic cross section.
double XS_inel(double Tp);
// 1 pion production cross section valid for Tp <= 2 GeV.
double XS_1pi(double Tp);
// 2 pion production cross section valid for Tp <= 2 GeV.
double XS_2pi(double Tp);
//=======================================================
// Multiplicities valid above Tp > 2 GeV, for:
// Geant4
double multip_pi0_Geant4(double Tp);
// Pythia8
double multip_pi0_Pythia8(double Tp);
// SIBYLL
double multip_pi0_SIBYLL(double Tp);
// QGSJET
double multip_pi0_QGSJET(double Tp);
//=======================================================
// Production cross section (=inelastic cross section*multiplicity) 
// valid above Tp > 2 GeV, for:
// Geant4
double XS_pi0_Geant4(double Tp);
// Pythia8
double XS_pi0_Pythia8(double Tp);
// SIBYLL
double XS_pi0_SIBYLL(double Tp);
// QGSJET
double XS_pi0_QGSJET(double Tp);
//=======================================================
// The peak value of the gamma-ray differential cross section in [mb/GeV], for:
// Geant4
double Amax_Geant4(double Tp);
// Pythia8
double Amax_Pythia8(double Tp);
// SIBYLL
double Amax_SIBYLL(double Tp);
// QGSJET
double Amax_QGSJET(double Tp);
//=======================================================
// The shape of the gamma-ray differential cross section function, for:
// Geant4
double F_Geant4(double Tp, double Egamma);
// Pythia8
double F_Pythia8(double Tp, double Egamma);
// SIBYLL
double F_SIBYLL(double Tp, double Egamma);
// QGSJET
double F_QGSJET(double Tp, double Egamma);
//=======================================================
// Gamma-ray differential cross section function in [mb/GeV], for:
// Geant4
double dXSdEg_Geant4(double Tp, double Egamma);
// Pythia8
double dXSdEg_Pythia8(double Tp, double Egamma);
// SIBYLL
double dXSdEg_SIBYLL(double Tp, double Egamma);
// QGSJET
double dXSdEg_QGSJET(double Tp, double Egamma);
//=======================================================
#endif
