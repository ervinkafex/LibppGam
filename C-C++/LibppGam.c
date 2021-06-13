/*
Writen by Ervin Kafexhiu, Oct 2014 and updated in May 2021
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


//=======================================================
// Tp - is the proton kinetic energy in [GeV].
// Egamma - is the gamma-ray energy in [GeV].
//=======================================================

#include <math.h>

#include "LibppGam.h"

//ooooooooooo0000000000oooooooooo0000000000ooooooooo#
//ooooooooooo0000000000oooooooooo0000000000ooooooooo#
//=======================================================
double Epi0_max_LAB(double Tp){
    //
	//This function calculates the maximum pi0 energy that
	//is allowed by the kinematics the LAB frame.
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]	
	//E_pi_max_LAB - is in [GeV]
    //
	double s, gamma_CM, Beta_CM; 
	double E_pi_CM, P_pi_CM;
	
	s = 2.0*m_p*(Tp + 2.0*m_p);
	gamma_CM = (Tp + 2.0*m_p)/sqrt(s);
	E_pi_CM = (s - 4.0*m_p*m_p + m_pi*m_pi)/(2.0*sqrt(s));
	P_pi_CM = sqrt(E_pi_CM*E_pi_CM - m_pi*m_pi);
	Beta_CM = sqrt(1.0 - 1.0/(gamma_CM*gamma_CM));
	return gamma_CM*(E_pi_CM + P_pi_CM*Beta_CM); // in GeV
}

double Egamma_max(double Tp){
    //
	//This function calculates the maximum gamma-ray energy
	//allowed by the kinematics in the LAB frame.
	//Tp and Egamma_max are in GeV.
	//
	double gamma_pi_LAB, Beta_pi_LAB;
	gamma_pi_LAB = Epi0_max_LAB(Tp)/m_pi;
	Beta_pi_LAB = sqrt(1.0 - 1.0/(gamma_pi_LAB*gamma_pi_LAB));
	return (m_pi/2.0)*gamma_pi_LAB*(1.0+Beta_pi_LAB); // in GeV
}


//=======================================================

double XS_inel(double Tp){
    //
	//This function calculates the pp total inelastic cross section
	//It originates from fitting the PDG data including TOTEM @ LHC.
    //
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//sigma_inel - inelastic cross section in [mb]
	//
	double LX, Threshold, ksi;
	LX  = log(Tp/Tp_th);
	ksi = 1.0 - pow( (Tp_th/Tp), 1.9);
	Threshold = (ksi > 0.0) ? ksi : 0.0;
	return (30.7 - 0.96*LX + 0.18*LX*LX)*(Threshold*Threshold*Threshold); //mb
}

double XS_1pi(double Tp){
    //
	//This function calculates the one pi0 production cross section.
	//It is valid for Tp <= 2 GeV. The channel included is pp->pp(pi0)
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//sigma_1_pi - is one pi0 production cross section in [mb]
	//
    	double M_RES     = 1.1883; // resonance effective mass in GeV
    	double Gamma_RES = 0.2264; // resonance effective width in GeV
    	double sigma_0   = 7.66e-3;// mb
    	double s, sqrts, X, eta, g_RES, K_RES, f_BW, XS1pi;

        if ((Tp <= Tp_th) && (Tp > 2.0)){XS1pi = 0.0;} // mb
        else {
            s = 2.0*m_p*(Tp + 2.0*m_p);
            sqrts = sqrt(s);
            X = sqrts - m_p;
            eta=sqrt((s-m_pi*m_pi-4.0*m_p*m_p)*(s-m_pi*m_pi-4.0*m_p*m_p) -16.0*m_pi*m_pi*m_p*m_p)/(2*m_pi*sqrts);
            //----------
            g_RES = sqrt((M_RES*M_RES)*(M_RES*M_RES + Gamma_RES*Gamma_RES));
            K_RES = sqrt(8.0)*M_RES*Gamma_RES*g_RES/(pi*sqrt(M_RES*M_RES + g_RES));
            f_BW = m_p*K_RES/((X*X-M_RES*M_RES)*(X*X-M_RES*M_RES)+ M_RES*M_RES*Gamma_RES*Gamma_RES);
	        //----------
            XS1pi = sigma_0*pow(eta,1.95)*(1.0+eta+pow(eta,5.0))*pow(f_BW,1.86); //mb
        }
        return XS1pi; // mb
}

double XS_2pi(double Tp){
    //
	//This function calculates the 2 pion production cross section.
	//It is valid Tp <= 2 GeV. The channels included here are
	//1) p+p -> p+n +(pi+)+(pi0)
	//2) p+p ->  D  +(pi+)+(pi0)
	//3) p+p -> p+p +2(pi0)
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//sigma_2_pi - is the sum of the pi0 cross sections from the above channels in [mb]
	//
	double XS2pi;
	if ((Tp < 0.56) && (Tp > 2.0)){XS2pi = 0.0;}
	else {
	    XS2pi =  5.7/(1.0 + exp(-9.3*(Tp-1.4))); //mb
	}      
    return XS2pi; //mb
}

//=======================================================
// Multiplicities valid above Tp > 2 GeV, for:
// Geant4
double multip_pi0_Geant4(double Tp){
    //
	//This function calculates the average pi0 production multiplicity given by Geant4.
	//This function is valid for TP>=1 GeV
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//multip_pi0 - is the average pi0 production multiplicity
	//
        double a_1,a_2,a_3,a_4,a_5;
        double Qp, xi_p, multip_pi0;
        if (Tp <= 2.0){multip_pi0 = 0.0;}  //GeV
    	else if ((Tp > 2.0) && (Tp < 5.0)){
    	    Qp = (Tp - Tp_th)/m_p;
            multip_pi0 = -6.0e-3 + 0.237*Qp - 0.023*Qp*Qp;
    	}else{
    	 // Geant4:
	     // --------------
			a_1 = 0.728;
			a_2 = 0.596;
			a_3 = 0.491;
			a_4 = 0.2503;
			a_5 = 0.117;
			xi_p = (Tp - 3.0)/m_p;
			multip_pi0 = a_1*pow(xi_p,a_4)*(1.0 + exp(-a_2*pow(xi_p,a_5)))*(1.0 - exp(-a_3*pow(xi_p,0.25)));
        // ---------------
    	}
	    return multip_pi0;
}
// Pythia8
double multip_pi0_Pythia8(double Tp){
    //
	//This function calculates the average pi0 production multiplicity given
	//by Geant4 for Tp <= 50 GeV and Pythia 8 for Tp > 50 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//multip_pi0 - is the average pi0 production multiplicity
	//
	    double a_1,a_2,a_3,a_4,a_5;
	    double xi_p, multip_pi0;
    	if (Tp <= 50.0){multip_pi0 = multip_pi0_Geant4(Tp);}
    	else {
    	 // Pythia8:
	    // ---------------
			a_1 =0.652 ;
			a_2 =0.0016;
			a_3 =0.488 ;
			a_4 =0.1928;
			a_5 =0.483 ;
	    
			xi_p = (Tp - 3.0)/m_p;
			multip_pi0 = a_1*pow(xi_p,a_4)*(1.0 + exp(-a_2*pow(xi_p,a_5)))*(1.0 - exp(-a_3*pow(xi_p,0.25)));
        // ---------------
    	}
        return multip_pi0;
}

// SIBYLL
double multip_pi0_SIBYLL(double Tp){
    //
	//This function calculates the average pi0 production multiplicity given
	//by Geant4 for Tp <= 100 GeV and SIBYLL for Tp > 100 GeV.
    //
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//multip_pi0 - is the average pi0 production multiplicity
	//
	    double a_1,a_2,a_3,a_4,a_5;
	    double xi_p, multip_pi0;
    	if (Tp <= 100.0){multip_pi0 = multip_pi0_Geant4(Tp);}
    	else {
    	 // SIBYLL:
	    // --------------
			a_1 = 5.436;
			a_2 = 0.254;
			a_3 = 0.072;
			a_4 = 0.075;
			a_5 = 0.166;
			xi_p = (Tp - 3.0)/m_p;
			multip_pi0 = a_1*pow(xi_p,a_4)*(1.0 + exp(-a_2*pow(xi_p,a_5)))*(1.0 - exp(-a_3*pow(xi_p,0.25)));
        // --------------
    	}
        return multip_pi0;
}

// QGSJET
double multip_pi0_QGSJET(double Tp){
    //
	//This function calculates the average pi0 production multiplicity given
	//by Geant4 for Tp <= 100 GeV and QGSJET for Tp > 100 GeV.
    //
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//multip_pi0 - is the average pi0 production multiplicity
	//
	    double a_1,a_2,a_3,a_4,a_5;
	    double xi_p, multip_pi0;
    	if (Tp <= 100.0){multip_pi0 = multip_pi0_Geant4(Tp);}
    	else {
    	 //QGSJET:
	    // ---------------
			a_1 = 0.908;
			a_2 = 0.0009;
			a_3 = 6.089;
			a_4 = 0.176;
			a_5 = 0.448;
			xi_p = (Tp - 3.0)/m_p;
			multip_pi0 = a_1*pow(xi_p,a_4)*(1.0 + exp(-a_2*pow(xi_p,a_5)))*(1.0 - exp(-a_3*pow(xi_p,0.25)));
        // ----------------
    	}
        return multip_pi0;
}

//=======================================================
// Production cross section (=inelastic cross section*multiplicity) 
// valid above Tp > 2 GeV, for:
// Geant4
double XS_pi0_Geant4(double Tp){
    //
	//This function calculates the pi0 production cross section
	//by using experimental data for Tp <= 2 GeV, and Geant4
	//multiplicity for Tp > 2 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//XS_pi0 - is the pi0 production cross section in [mb]
	//
    return XS_1pi(Tp) + XS_2pi(Tp) + XS_inel(Tp)*multip_pi0_Geant4(Tp);  
}

// Pythia8
double XS_pi0_Pythia8(double Tp){
    //
	//This function calculates the pi0 production cross section
	//by using experimental data for Tp <= 2 GeV, Geant4
	//multiplicity for 2 < Tp <= 50 GeV and Pythia8 multiplicity for Tp > 50 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//XS_pi0 - is the pi0 production cross section in [mb]
	//
    return XS_1pi(Tp) + XS_2pi(Tp) + XS_inel(Tp)*multip_pi0_Pythia8(Tp);
}

// SIBYLL
double XS_pi0_SIBYLL(double Tp){
    //
	//This function calculates the pi0 production cross section
	//by using experimental data for Tp <= 2 GeV, Geant4 multiplicity
	//for 2 < Tp <= 100 GeV and SIBYLL multiplicity for Tp > 100 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//XS_pi0 - is the pi0 production cross section in [mb]
	//
    return XS_1pi(Tp) + XS_2pi(Tp) + XS_inel(Tp)*multip_pi0_SIBYLL(Tp);
}

// QGSJET
double XS_pi0_QGSJET(double Tp){
    //
	//This function calculates the pi0 production cross section
	//by using experimental data for Tp <= 2 GeV, Geant4 multiplicity
	//for 2 < Tp <= 100 GeV and QGSJET multiplicity for Tp > 100 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//XS_pi0 - is the pi0 production cross section in [mb]
	//
    return XS_1pi(Tp) + XS_2pi(Tp) + XS_inel(Tp)*multip_pi0_QGSJET(Tp);
}

//=======================================================
// The peak value of the gamma-ray differential cross section in [mb/GeV], for:
// Geant4
double Amax_Geant4(double Tp){
    //
	//This function calculates the peak value of the gamma-ray
	//differential cross section. For Tp < 1 GeV is the fit from
	//the experimental data, and for Tp >= 1 GeV is based on Geant4.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	//
	
	double theta_p = Tp/m_p;
	double Ltheta_p= log(theta_p);
	double Amax;
	double b_1,b_2,b_3;

	if (Tp <= Tp_th) {Amax = 0.0;}
	else if ((Tp > Tp_th) && (Tp < 1.0)){Amax = 5.9*XS_pi0_Geant4(Tp)/Epi0_max_LAB(Tp);}
	else if ((Tp >= 1.0)  && (Tp < 5.0)){
	    b_1 =  9.53;
		b_2 = -0.52;
		b_3 =  0.054;
		Amax = b_1*pow(theta_p,b_2)*exp(b_3*Ltheta_p*Ltheta_p)*XS_pi0_Geant4(Tp)/m_p;
	}
	else{
	    b_1 =  9.13;
		b_2 = -0.35;
		b_3 =  9.7e-3;
		Amax = b_1*pow(theta_p,b_2)*exp(b_3*Ltheta_p*Ltheta_p)*XS_pi0_Geant4(Tp)/m_p;
	}
	return Amax;
}

// Pythia8
double Amax_Pythia8(double Tp){
	//
	//This function calculates the peak value of the gamma-ray
	//differential cross section. For Tp < 1 GeV is the fit from
	//the experimental data, for 1 <= Tp <= 50 GeV is based on Geant4
	//and for Tp > 50 GeV is based on Pythia8.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	//
	
    double theta_p = Tp/m_p;
	double Ltheta_p= log(theta_p);
    double Amax;
    double b_1,b_2,b_3;
    
	if(Tp <= 50.0) {Amax = Amax_Geant4(Tp);}
	else{
	    b_1 =  9.06;
		b_2 = -0.3795;
		b_3 =  0.01105;
		Amax = b_1*pow(theta_p,b_2)*exp(b_3*Ltheta_p*Ltheta_p)*XS_pi0_Pythia8(Tp)/m_p;
	}	
	return Amax;
}

// SIBYLL
double Amax_SIBYLL(double Tp){
    //
	//This function calculates the peak value of the gamma-ray
	//differential cross section. For Tp < 1 GeV is the fit from
	//the experimental data, for 1 <= Tp <= 100 GeV is based on Geant4
	//and for Tp > 100 GeV is based on SIBYLL.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	//
	
    double theta_p = Tp/m_p;
	double Ltheta_p= log(theta_p);
	double Amax;
    double b_1,b_2,b_3;

	if(Tp <= 100.0) {Amax = Amax_Geant4(Tp);}
	else{
	    b_1 =  10.77;
		b_2 = -0.412;
		b_3 =  0.01264;
		Amax = b_1*pow(theta_p,b_2)*exp(b_3*Ltheta_p*Ltheta_p)*XS_pi0_SIBYLL(Tp)/m_p;
	}	
	return Amax;
}
// QGSJET
double Amax_QGSJET(double Tp){
    //
	//This function calculates the peak value of the gamma-ray
	//differential cross section. For Tp < 1 GeV is the fit from
	//the experimental data, for 1 <= Tp <= 100 GeV is based on Geant4
	//and for Tp > 100 GeV is based on QGSJET.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	//
	
    double theta_p = Tp/m_p;
	double Ltheta_p= log(theta_p);
	double Amax;
    double b_1,b_2,b_3;

	if(Tp <= 100.0) {Amax = Amax_Geant4(Tp);}
	else{
	    b_1 =  13.16;
		b_2 = -0.4419;
		b_3 =  0.01439;
		Amax = b_1*pow(theta_p,b_2)*exp(b_3*Ltheta_p*Ltheta_p)*XS_pi0_QGSJET(Tp)/m_p;
	}	
		
	return Amax;
}

//=======================================================
// The shape of the gamma-ray differential cross section function, for:
// Geant4
double F_Geant4(double Tp, double Egamma){
    //
	//This function calculates the shape of the gamma-ray
	//differential cross section function for a specific Tp.
	//This function includes the experimental data 
	//for Tp < 1 GeV and Geant4 for Tp >= 1 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//FF - is the shape of the gamma-ray spectrum, it is unitless.
	//
	
	double Y  = Egamma + m_pi*m_pi/(4.0*Egamma);
	double Egmax = Egamma_max(Tp);
	double Y0 = Egmax + m_pi*m_pi/(4.0*Egmax);
	double X  = (Y-m_pi)/(Y0-m_pi);
	double theta = Tp/m_p;
	double kappa = 3.29 - 0.2*pow(theta, -1.5);
	// --------------
	double q  = (Tp -1.0)/m_p;
	double C  = 3.0*m_pi/Y0;
	double mu, beta, gamma, FF;
	
	if ((X >= 0.0) && (X < 1.0)){
	    if((Tp > Tp_th) && (Tp<1.0)){FF = pow(1.0-X, kappa);}
		else if ((Tp>=1.0) && (Tp<=4.0)){
		    mu = 1.25*pow(q, 1.25) * exp(-1.25*q);
			beta = mu + 2.45;
			gamma= mu + 1.45;
			FF = pow(1.0-X, beta)/pow(1.0+X/C, gamma);
		}
		else if ((Tp>4.0) && (Tp<=20.0)){
		    mu = 1.25*pow(q, 1.25) * exp(-1.25*q);
			beta = 1.5*mu + 4.95;
			gamma= mu + 1.5;
			FF = pow(1.0-X, beta)/pow(1.0+X/C, gamma);
		}
		else if ((Tp>20.0) && (Tp<=100.0)){
		    FF = pow(1-sqrt(X), 4.2)/(1.0+X/C);
		}
		else if (Tp>100.0){
		    FF = pow(1-sqrt(X), 4.9)/(1.0+X/C);
		}
		else{FF = 0.0;}
	}
	else{FF = 0.0;}	
	return FF;
}

// Pythia8
double F_Pythia8(double Tp, double Egamma){
    //
	//This function calculates the shape of the gamma-ray
	//differential cross section function for a specific Tp.
	//This function includes the experimental data for
	//Tp < 1 GeV Geant4 for 1 <= Tp <= 50 GeV and Pythia8 for Tp > 50 GeV
	//	
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//FF - is the shape of the gamma-ray spectrum, it is unitless.
	//
	
	double Y  = Egamma + m_pi*m_pi/(4.0*Egamma);
	double Egmax = Egamma_max(Tp);
	double Y0 = Egmax + m_pi*m_pi/(4.0*Egmax);
	double X  = (Y-m_pi)/(Y0-m_pi);
	double C  = 3.5*m_pi/Y0;
	double FF;
	
	if ((X >= 0.0) && (X < 1.0) && (Tp > 50)){
	    FF = pow(1-sqrt(X), 4.0)/(1.0+X/C);
	} 	
	else{FF = F_Geant4(Tp, Egamma);}
	return FF;
}

// SIBYLL
double F_SIBYLL(double Tp, double Egamma){
    //
	//This function calculates the shape of the gamma-ray
	//differential cross section function for a specific Tp.
	//This function includes the experimental data for
	//Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and SIBYLL for Tp > 100 GeV
	//	
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//FF - is the shape of the gamma-ray spectrum, it is unitless.
	//
	
	double Y  = Egamma + m_pi*m_pi/(4.0*Egamma);
	double Egmax = Egamma_max(Tp);
	double Y0 = Egmax + m_pi*m_pi/(4.0*Egmax);
	double X  = (Y-m_pi)/(Y0-m_pi);
	double C  = 3.55*m_pi/Y0;
	double FF;
	
	if ((X >= 0.0) && (X < 1.0) && (Tp > 100)){
	    FF = pow(1-sqrt(X), 3.6)/(1.0+X/C);
	} 	
	else{FF = F_Geant4(Tp, Egamma);}
	return FF;
}

// QGSJET
double F_QGSJET(double Tp, double Egamma){
    //
	//This function calculates the shape of the gamma-ray
	//differential cross section function for a specific Tp.
	//This function includes the experimental data for
	//Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and QGSJET for Tp > 100 GeV
	//	
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//FF - is the shape of the gamma-ray spectrum, it is unitless.
	//
	
	double Y  = Egamma + m_pi*m_pi/(4.0*Egamma);
	double Egmax = Egamma_max(Tp);
	double Y0 = Egmax + m_pi*m_pi/(4.0*Egmax);
	double X  = (Y-m_pi)/(Y0-m_pi);
	double C  = 3.55*m_pi/Y0;
	double FF;
	
	if ((X >= 0.0) && (X < 1.0) && (Tp > 100)){
	    FF = pow(1-sqrt(X), 4.5)/(1.0+X/C);
	} 	
	else{FF = F_Geant4(Tp, Egamma);}
	return FF;
}

//=======================================================
// Gamma-ray differential cross section function in [mb/GeV], for:
// Geant4
double dXSdEg_Geant4(double Tp, double Egamma){
    //
	//This function calculates the pp->pi0->gamma-ray
	//differential cross section function.
	//This function includes the experimental data 
	//for Tp < 1 GeV and Geant4 for Tp >= 1 GeV.
	//
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//The value that is returned is in [mb/GeV]
	// 
	return Amax_Geant4(Tp)*F_Geant4(Tp, Egamma);
}

// Pythia8
double dXSdEg_Pythia8(double Tp, double Egamma){
    //
	//This function calculates the pp->pi0->gamma-ray
	//differential cross section function.
	//This function includes the experimental data for
	//Tp < 1 GeV Geant4 for 1 <= Tp <= 50 GeV and Pythia8 for Tp > 50 GeV
	//	
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//The value that is returned is in [mb/GeV]
	//
	return Amax_Pythia8(Tp)*F_Pythia8(Tp, Egamma);
}

// SIBYLL
double dXSdEg_SIBYLL(double Tp, double Egamma){
    //
	//This function calculates the pp->pi0->gamma-ray
	//differential cross section function.
	//This function includes the experimental data for
	//Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and SIBYLL for Tp > 100 GeV
	//	
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//The value that is returned is in [mb/GeV]
	//
	return Amax_SIBYLL(Tp)*F_SIBYLL(Tp, Egamma);
}

// QGSJET
double dXSdEg_QGSJET(double Tp, double Egamma){
    //
	//This function calculates the pp->pi0->gamma-ray
	//differential cross section function.
	//This function includes the experimental data for
	//Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and QGSJET for Tp > 100 GeV
	//	
	//Tp - is the proton kinetic energy in the LAB frame in [GeV]
	//The value that is returned is in [mb/GeV]
	//
	return Amax_QGSJET(Tp)*F_QGSJET(Tp, Egamma);
}

//=======================================================
