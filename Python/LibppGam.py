# -*- coding: utf-8 -*-

""" 
Writen by Pol Bordas & Ervin Kafexhiu, Oct 2014.
Updated by Ervin Kafexhiu, May 2021 (compliant with python 2&3 versions).
===========================================================
     _     _ _                 _____
    | |   (_) |               |  __ \
    | |    _| |__  _ __  _ __ | |  \/ __ _ _ __ ___
    | |   | | '_ \| '_ \| '_ \| | __ / _` | '_ ` _ \
    | |___| | |_) | |_) | |_) | |_\ \ (_| | | | | | |
    \_____/_|_.__/| .__/| .__/ \____/\__,_|_| |_| |_|
                  | |   | |
                  |_|   |_|
===========================================================
If used please reference to Phys.Rev. D90 (2014) 12, 123014 
(astro-ph/1406.7369), "Parametrization of gamma-ray 
production cross-sections for pp interactions in a broad 
proton energy range from the kinematic threshold to PeV energies" 
by Ervin Kafexhiu, Felix Aharonian, Andrew M. Taylor, Gabriela S. Vila
===========================================================
This library contains the parametrization for different 
quantities for the pp->pi0->2gamma.
The basic units in each functions are GeV and mb.
"""

#ooooooooooo0000000000oooooooooo0000000000ooooooooo#
#ooooooooooo0000000000oooooooooo0000000000ooooooooo#

# import math library
import math

# definition of some constants

m_p  = 0.938272 # GeV, proton mass (taken from PDG)
m_pi = 0.134976 # GeV, pi0 mass (taken from PDG)

# the proton threshold energy in the LAB frame
Tp_th = 2.0*m_pi + m_pi*m_pi/(2.0*m_p) # in GeV

#ooooooooooo0000000000oooooooooo0000000000ooooooooo#
#ooooooooooo0000000000oooooooooo0000000000ooooooooo#

############################
def Epi0_max_LAB(Tp):
    """
	This function calculates the maximum pi0 energy that
	is allowed by the kinematics the LAB frame.
	Tp - is the proton kinetic energy in the LAB frame in [GeV]	
	Epi_maxLAB - is in [GeV]
	"""
    s = 2.0*m_p*(Tp + 2.0*m_p)
    gamma_CM = (Tp + 2.0*m_p)/(s**0.5)
    E_pi_CM  = (s - 4.0*m_p**2.0 + m_pi**2.0)/(2.0*s**0.5)
    P_pi_CM  = ( E_pi_CM**2.0 - m_pi**2.0 )**0.5
    Beta_CM  = ( 1.0 - gamma_CM**(-2.0) )**0.5
    Epi_maxLAB = gamma_CM*(E_pi_CM + P_pi_CM*Beta_CM) # in GeV
    return Epi_maxLAB

############################
def Egamma_max(Tp):
    """
	This function calculates the maximum gamma-ray energy
	allowed by the kinematics in the LAB frame.
	Tp and Eg_max are in GeV.
	"""
    gamma_pi_LAB = Epi0_max_LAB(Tp)/m_pi
    Beta_pi_LAB = (1.0 - gamma_pi_LAB**(-2.0))**0.5
    #Eg_min=(m_pi/2.0)*gamma_pi_LAB*(1.0-Beta_pi_LAB)
    Eg_max = (m_pi/2.0)*gamma_pi_LAB*(1.0+Beta_pi_LAB)
    return Eg_max

############################
def sigma_inel(Tp):
    """
	This function calculates the pp total inelastic cross section
	It originates from fitting the PDG data including TOTEM @ LHC.

	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS - inelastic cross section in [mb]
	"""
    LX = math.log(Tp/Tp_th)
    Threshold = max( 0.0, (1.0-(Tp_th/Tp)**1.9) )
    XS = (30.7 - 0.96*LX+ 0.18*LX**2.0)*(Threshold**3.0)
    return XS
		
############################
def sigma_1pi(Tp):
    """
	This function calculates the one pi0 production cross section.
	It is valid for Tp <= 2 GeV. The channel included is pp->pp(pi0)
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS_1_pi - is one pi0 production cross section in [mb]
	"""
    M_RES     = 1.1883 # resonance effective mass in GeV
    Gamma_RES = 0.2264 # resonance effective width in GeV
    sigma_0   = 7.66e-3# mb
    if (Tp > Tp_th) and (Tp <= 2.0):
        s = 2.0*m_p*(Tp + 2.0*m_p)
        X = s**0.5 - m_p
        eta=(( (s-m_pi**2.0 -4.0*m_p**2.0)**2.0 -16.0*(m_pi**2.0)*(m_p**2.0))**0.5)/(2*m_pi*(s**0.5))
        #----------
        g_RES = ((M_RES**2.0)*(M_RES**2.0 + Gamma_RES**2.0))**0.5
        K_RES = (8.0**0.5)*M_RES*Gamma_RES*g_RES/(math.pi*(M_RES**2.0 + g_RES)**0.5)
        f_BW = m_p*K_RES/((X**2.0 - M_RES**2.0)**2.0 +(M_RES*Gamma_RES)**2.0)
        #--------
        XS_1_pi = sigma_0*(eta**1.95)*(1.0 +eta +eta**5.0)*(f_BW**1.86) #mb
    else:
        XS_1_pi =  0.0 # mb
    return XS_1_pi

############################
def sigma_2pi(Tp):
    """
	This function calculates the 2 pion production cross section.
	It is valid Tp <= 2 GeV. The channels included here are
	1) p+p -> p+n +(pi+)+(pi0)
	2) p+p ->  D  +(pi+)+(pi0)
	3) p+p -> p+p +2(pi0)
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS_2_pi - is the sum of the pi0 cross sections from the above channels in [mb]
	"""
    if (Tp >= 0.56) and (Tp <= 2.0):
        XS_2_pi = 5.7/(1.0 + math.exp(-9.3*(Tp-1.4))) #mb
    else:
        XS_2_pi = 0.0
    return XS_2_pi #mb

############################
def multip_pi0_Geant4(Tp):
    """
	This function calculates the average pi0 production multiplicity given by Geant4.
	This function is valid for TP>=1 GeV
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	multip_pi0 - is the average pi0 production multiplicity
	"""
    if Tp <= 2.0:
        multip_pi0 = 0.0
    elif (Tp > 2.0) and (Tp < 5.0):
        Qp = (Tp - Tp_th)/m_p
        multip_pi0 = -6.0e-3 + 0.237*Qp - 0.023*(Qp**2.0)
    else:
        # Geant4:
        # -----------
        a_1 = 0.728
        a_2 = 0.596
        a_3 = 0.491
        a_4 = 0.2503
        a_5 = 0.117
        # -----------
        xi_p = (Tp - 3.0)/m_p
        multip_pi0 = a_1*(xi_p**a_4)*(1.0 + math.exp(-a_2*(xi_p**a_5)))*(1.0 - math.exp(-a_3*(xi_p**0.25)))
    return multip_pi0
        
############################
def multip_pi0_Pythia8(Tp):
    """
	This function calculates the average pi0 production multiplicity given
	by Geant4 for Tp <= 50 GeV and Pythia 8 for Tp > 50 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	multip_pi0 - is the average pi0 production multiplicity
	"""
    if Tp <= 50.0: #GeV
        # use Geant4:
        multip_pi0 = multip_pi0_Geant4(Tp)
    else:
        # Pythia8:
        # -----------
        a_1 =0.652 
        a_2 =0.0016
        a_3 =0.488 
        a_4 =0.1928
        a_5 =0.483
        # -----------
        xi_p = (Tp - 3.0)/m_p
        multip_pi0 = a_1*(xi_p**a_4)*(1.0 + math.exp(-a_2*(xi_p**a_5)))*(1.0 - math.exp(-a_3*(xi_p**0.25)))
    return multip_pi0 

############################
def multip_pi0_SIBYLL(Tp):
    """
	This function calculates the average pi0 production multiplicity given
	by Geant4 for Tp <= 100 GeV and SIBYLL for Tp > 100 GeV.

	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	multip_pi0 - is the average pi0 production multiplicity
	"""
    if Tp <= 100.0: #GeV
        # use Geant4:
        multip_pi0 = multip_pi0_Geant4(Tp)
    else:
        # SIBYLL:
        ## -----------
        a_1 = 5.436
        a_2 = 0.254
        a_3 = 0.072
        a_4 = 0.075
        a_5 = 0.166
	    # -----------
        xi_p = (Tp - 3.0)/m_p
        multip_pi0 = a_1*(xi_p**a_4)*(1.0 + math.exp(-a_2*(xi_p**a_5)))*(1.0 - math.exp(-a_3*(xi_p**0.25)))
    return multip_pi0 

############################
def multip_pi0_QGSJET(Tp):
    """
	This function calculates the average pi0 production multiplicity given
	by Geant4 for Tp <= 100 GeV and QGSJET for Tp > 100 GeV.

	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	multip_pi0 - is the average pi0 production multiplicity
	"""
    if Tp <= 100.0:
        # use Geant4:
        multip_pi0 = multip_pi0_Geant4(Tp)
    else:
        # QGSJET:
        # -----------
        a_1 = 0.908
        a_2 = 0.0009
        a_3 = 6.089
        a_4 = 0.176
        a_5 = 0.448
        # -----------
        xi_p = (Tp - 3.0)/m_p
        multip_pi0 = a_1*(xi_p**a_4)*(1.0 + math.exp(-a_2*(xi_p**a_5)))*(1.0 - math.exp(-a_3*(xi_p**0.25)))
    return multip_pi0

############################
def sigma_pi_Geant4(Tp):
    """
	This function calculates the pi0 production cross section
	by using experimental data for Tp <= 2 GeV, and Geant4
	multiplicity for Tp > 2 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS_pi0 - is the pi0 production cross section in [mb]
	"""
    return sigma_1_pi(Tp) + sigma_2_pi(Tp) + sigma_inel(Tp)*multip_pi0_Geant4(Tp)

############################ 
def sigma_pi_Pythia8(Tp):
    """
	This function calculates the pi0 production cross section
	by using experimental data for Tp <= 2 GeV, Geant4
	multiplicity for 2 < Tp <= 50 GeV and Pythia8 multiplicity for Tp > 50 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS_pi0 - is the pi0 production cross section in [mb]
	"""
    return sigma_1_pi(Tp) + sigma_2_pi(Tp) + sigma_inel(Tp)*multip_pi0_Pythia8(Tp)

############################ 
def sigma_pi_SIBYLL(Tp):
    """
	This function calculates the pi0 production cross section
	by using experimental data for Tp <= 2 GeV, Geant4 multiplicity
	for 2 < Tp <= 100 GeV and SIBYLL multiplicity for Tp > 100 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS_pi0 - is the pi0 production cross section in [mb]
	"""
    return sigma_1_pi(Tp) + sigma_2_pi(Tp) + sigma_inel(Tp)*multip_pi0_SIBYLL(Tp)
    
############################ 
def sigma_pi_QGSJET(Tp):
    """
	This function calculates the pi0 production cross section
	by using experimental data for Tp <= 2 GeV, Geant4 multiplicity
	for 2 < Tp <= 100 GeV and QGSJET multiplicity for Tp > 100 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	XS_pi0 - is the pi0 production cross section in [mb]
	"""
    return sigma_1_pi(Tp) + sigma_2_pi(Tp) + sigma_inel(Tp)*multip_pi0_QGSJET(Tp)

############################
def Amax_Geant4(Tp):
    """
	This function calculates the peak value of the gamma-ray
	differential cross section. For Tp < 1 GeV is the fit from
	the experimental data, and for Tp >= 1 GeV is based on Geant4.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	"""
    theta_p = Tp/m_p
    Ltheta_p= math.log(theta_p)
    
    if Tp <= Tp_th:
        Amax = 0.0
    elif (Tp > Tp_th) and (Tp < 1.0):
        Amax = 5.9*sigma_pi_Geant4(Tp)/E_pi_max_LAB(Tp)
    elif (Tp >= 1.0) and (Tp < 5.0):
        b_1 =  9.53
        b_2 = -0.52
        b_3 =  0.054
        Amax = b_1*(theta_p**b_2)*math.exp(b_3*Ltheta_p**2.0)*sigma_pi_Geant4(Tp)/m_p
    else:
        b_1 =  9.13
        b_2 = -0.35
        b_3 =  9.7e-3
        Amax = b_1*(theta_p**b_2)*math.exp(b_3*Ltheta_p**2.0)*sigma_pi_Geant4(Tp)/m_p
    return Amax

############################
def Amax_Pythia8(Tp):
    """
	This function calculates the peak value of the gamma-ray
	differential cross section. For Tp < 1 GeV is the fit from
	the experimental data, for 1 <= Tp <= 50 GeV is based on Geant4
	and for Tp > 50 GeV is based on Pythia8.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	"""
    theta_p = Tp/m_p
    Ltheta_p= math.log(theta_p)
    if Tp <= 50.0:
        Amax = Amax_Geant4(Tp)
    else:
        b_1 =  9.06
        b_2 = -0.3795
        b_3 =  0.01105
        Amax = b_1*(theta_p**b_2)*math.exp(b_3*Ltheta_p**2.0)*sigma_pi_Pythia8(Tp)/m_p
    return Amax

############################
def Amax_SIBYLL(Tp):
    """
	This function calculates the peak value of the gamma-ray
	differential cross section. For Tp < 1 GeV is the fit from
	the experimental data, for 1 <= Tp <= 100 GeV is based on Geant4
	and for Tp > 100 GeV is based on SIBYLL.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	"""
    theta_p = Tp/m_p
    Ltheta_p= math.log(theta_p)
    if Tp <= 100.0:
        Amax = Amax_Geant4(Tp)
    else:
        b_1 =  10.77
        b_2 = -0.412
        b_3 =  0.01264
        Amax = b_1*(theta_p**b_2)*math.exp(b_3*Ltheta_p**2.0)*sigma_pi_SIBYLL(Tp)/m_p
    return Amax

############################
def Amax_QGSJET(Tp):
    """
	This function calculates the peak value of the gamma-ray
	differential cross section. For Tp < 1 GeV is the fit from
	the experimental data, for 1 <= Tp <= 100 GeV is based on Geant4
	and for Tp > 100 GeV is based on QGSJET.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	Amax - is the peak of the gamma-ray differential cross section in [mb/GeV]
	"""
    theta_p = Tp/m_p
    Ltheta_p= math.log(theta_p)
    if Tp <= 100.0:
        Amax = Amax_Geant4(Tp)
    else:
        b_1 =  13.16
        b_2 = -0.4419
        b_3 =  0.01439
        Amax = b_1*(theta_p**b_2)*math.exp(b_3*Ltheta_p**2.0)*sigma_pi_QGSJET(Tp)/m_p
    return Amax

############################
def F_Geant4(Tp, Egamma):
    """
	This function calculates the shape of the gamma-ray
	differential cross section function for a specific Tp.
	This function includes the experimental data 
	for Tp < 1 GeV and Geant4 for Tp >= 1 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	FF - is the shape of the gamma-ray spectrum, it is unitless.
	"""
    Y  = Egamma + (m_pi**2.0)/(4.0*Egamma)
    Y0 = Egamma_max(Tp) + (m_pi**2.0)/(4.0*Egamma_max(Tp))
    X  = (Y-m_pi)/(Y0-m_pi)
    theta = Tp/m_p
    kappa = 3.29 - 0.2*theta**(-1.5)
    # --------------
    q  = (Tp -1.0)/m_p
    C  = 3.0*m_pi/Y0;
	# --------------
    if X >= 0.0 and X < 1.0:
        if (Tp > Tp_th) and (Tp<1.0):
            FF = (1.0-X)**kappa
        elif (Tp>=1.0) and (Tp<=4.0):
            mu = 1.25*(q**(1.25))*math.exp(-1.25*q)
            beta = mu + 2.45
            gamma= mu + 1.45
            FF = ((1.0-X)**beta)/((1.0+X/C)**gamma)
        elif Tp>4.0 and Tp<=20.0:
            mu = 1.25*(q**(1.25))*math.exp(-1.25*q)
            beta = 1.5*mu + 4.95
            gamma= mu + 1.5
            FF = ((1.0-X)**beta)/((1+X/C)**gamma)
        elif (Tp>20.0) and (Tp<=100.0):
            FF = ((1.0-X**0.5)**4.2)/(1.0+X/C)
        elif Tp>100.0:
            FF = ((1.0-X**0.5)**4.9)/(1.0+X/C)
        else:
            FF = 0.0
    else:
        FF = 0.0
    return FF

############################
def F_Pythia8(Tp, Egamma):
    """
	This function calculates the shape of the gamma-ray
	differential cross section function for a specific Tp.
	This function includes the experimental data for
	Tp < 1 GeV Geant4 for 1 <= Tp <= 50 GeV and Pythia8 for Tp > 50 GeV
		
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	FF - is the shape of the gamma-ray spectrum, it is unitless.
	"""
    Y  = Egamma + (m_pi**2.0)/(4.0*Egamma)
    Y0 = Egamma_max(Tp) + (m_pi**2.0)/(4.0*Egamma_max(Tp))
    X  = (Y-m_pi)/(Y0-m_pi)
    C  = 3.5*m_pi/Y0
    if (X >= 0.0) and (X < 1.0) and (Tp > 50):
        FF = ((1.0-X**0.5)**4.0)/(1+X/C)
    else:
        FF = F_Geant4(Tp, Egamma)
    return FF

############################
def F_SIBYLL(Tp, Egamma):
    """
	This function calculates the shape of the gamma-ray
	differential cross section function for a specific Tp.
	This function includes the experimental data for
	Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and SIBYLL for Tp > 100 GeV
		
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	FF - is the shape of the gamma-ray spectrum, it is unitless.
	"""
    Y  = Egamma + (m_pi**2.0)/(4.0*Egamma)
    Y0 = Egamma_max(Tp) + (m_pi**2.0)/(4.0*Egamma_max(Tp))
    X  = (Y-m_pi)/(Y0-m_pi)
    C  = 3.55*m_pi/Y0
    if (X >= 0.0) and (X < 1.0) and (Tp > 100):
        FF = ((1.0-X**0.5)**3.6)/(1.0+X/C)
    else:
        FF = F_Geant4(Tp, Egamma)
    return FF
 
############################ 
def F_QGSJET(Tp, Egamma):
    """
	This function calculates the shape of the gamma-ray
	differential cross section function for a specific Tp.
	This function includes the experimental data for
	Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and QGSJET for Tp > 100 GeV
		
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	FF - is the shape of the gamma-ray spectrum, it is unitless.
	"""
    Y  = Egamma + (m_pi**2.0)/(4.0*Egamma)
    Y0 = Egamma_max(Tp) + (m_pi**2.0)/(4.0*Egamma_max(Tp))
    X  = (Y-m_pi)/(Y0-m_pi)
    C  = 3.55*m_pi/Y0
    if (X >= 0.0) and (X < 1.0) and (Tp > 100):
        FF = ((1.0-X**0.5)**4.5)/(1.0+X/C)
    else:
        FF = F_Geant4(Tp, Egamma)
    return FF

############################
def dsigma_dEgamma_Geant4(Tp, Egamma):
    """
	This function calculates the pp->pi0->gamma-ray
	differential cross section function.
	This function includes the experimental data 
	for Tp < 1 GeV and Geant4 for Tp >= 1 GeV.
	
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	The value that is returned is in [mb/GeV]
	"""  
    return Amax_Geant4(Tp)*F_Geant4(Tp, Egamma)

############################
def dsigma_dEgamma_Pythia8(Tp, Egamma):
    """
	This function calculates the pp->pi0->gamma-ray
	differential cross section function.
	This function includes the experimental data for
	Tp < 1 GeV Geant4 for 1 <= Tp <= 50 GeV and Pythia8 for Tp > 50 GeV
		
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	The value that is returned is in [mb/GeV]
	"""
    return Amax_Pythia8(Tp)*F_Pythia8(Tp, Egamma)

############################
def dsigma_dEgamma_SIBYLL(Tp, Egamma):
    """
	This function calculates the pp->pi0->gamma-ray
	differential cross section function.
	This function includes the experimental data for
	Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and SIBYLL for Tp > 100 GeV
		
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	The value that is returned is in [mb/GeV]
	"""
    return Amax_SIBYLL(Tp)*F_SIBYLL(Tp, Egamma)

############################
def dsigma_dEgamma_QGSJET(Tp, Egamma):
    """
	This function calculates the pp->pi0->gamma-ray
	differential cross section function.
	This function includes the experimental data for
	Tp < 1 GeV Geant4 for 1 <= Tp <= 100 GeV and QGSJET for Tp > 100 GeV
		
	Tp - is the proton kinetic energy in the LAB frame in [GeV]
	The value that is returned is in [mb/GeV]
	"""
    return Amax_QGSJET(Tp)*F_QGSJET(Tp, Egamma)
############################
