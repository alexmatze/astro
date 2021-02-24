# --------------------------------------------------------
# The calculus can be looked up in the paper by D. W. Hogg
# https://arxiv.org/pdf/astro-ph/9905116.pdf
# --------------------------------------------------------

import numpy as np
from scipy.integrate import quad


q1 = input ("Are your values defined in the script? (Y/N - default: Yes): ") or "Y"


if str.lower(q1) == "yes" or str.lower(q1) == "y" or str.lower(q1) == "ye" or str.lower(q1) == "j" or str.lower(q1) == "ja":
    #define worldmodel:
    H0=71.0             # Hubble constant in km/s * Mpc
    omega_M = 0.27      # Matter-Term
    omega_lambda = 0.73 # Vacuum-Term
    omega_k = 0.0       # curvature

    #define redshift:
    z = 4.3             # cosmologic reshift
else:
    H0 =  input ("Hubble constant in km/s * Mpc? (default: 71.0): ") or "71.0"   # Hubble constant in km/s * Mpc
    omega_M =  input ("Omega_m? (default: 0.27): ") or "0.27"                    # Matter-Term
    omega_lambda =  input ("Omega_lambda? (default: 0.73): ") or "0.73"          # Vacuum-Term
    omega_k = input ("Omega_k? (default: 0.00): ") or "0.00"                     # curvature
    z = input ("cosmologic redshift? (default: 4.3): ") or "4.3"                 # cosmologic reshift
    H0=float(H0)
    omega_M=float(omega_M)
    omega_lambda=float(omega_lambda)
    omega_k=float(omega_k)
    z=float(z)



#transverse distance between 2 signals
def D_M(omega_k,omega_M,omega_Lambda,z,H0):
    if omega_k+omega_M+omega_Lambda != 1:
        print(r'make sure that $\Omega _k+\Omega _M + \Omega _{\Lambda}$ = 1')
        return 0
    H_0 = (H0/3.086)*10**(-19)
    D_h = 299792458/(H_0)
    E = lambda redshift,o_k,o_M,o_Lambda:1/np.sqrt(o_M*(1+redshift)**3+o_k*(1+redshift)**2+o_Lambda)
    D_c = D_h * quad(E,0,z,args=(omega_k,omega_M,omega_Lambda))[0]
    if omega_k > 0:
        Dc = D_h * (1/np.sqrt(omega_k)) * np.sinh(np.sqrt(omega_k)*D_c/D_h)
    if omega_k == 0:
        Dc = D_c
    if omega_k < 0:
        Dc = D_h * (1/np.sqrt(omega_k)) * np.sin(np.sqrt(omega_k)*D_c/D_h)
    return(Dc) # in m

#angular diameter distance
def D_A(omega_k,omega_M,omega_Lambda,z,H0):
    return D_M(omega_k, omega_M, omega_Lambda, z, H0)/(1+z) #returns m/rad

#distance line-of-sight
def D_c(omega_k,omega_M,omega_Lambda,z,H0):
    if omega_k+omega_M+omega_Lambda == 1:
        H_0 = (H0/3.086)*10**(-19)
        D_h = 299792458/(H_0)
        E = lambda redshift,o_k,o_M,o_Lambda:1/np.sqrt(o_M*(1+redshift)**3+o_k*(1+redshift)**2+o_Lambda)
        Dc = D_h * quad(E,0,z,args=(omega_k,omega_M,omega_Lambda))[0]
        return(Dc) # in m
    else:
        print(r'make sure that $\Omega _k+\Omega _M + \Omega _{\Lambda}$ = 1')
        return 0


#lookback time
def t_l(omega_k,omega_M,omega_Lambda,z,H0):
    if omega_k+omega_M+omega_Lambda != 1:
        print(r'make sure that $\Omega _k+\Omega _M + \Omega _{\Lambda}$ = 1')
        return 0
    H_0 = (H0/3.086)*10**(-19)
    t_h = 1/(H_0)
    E = lambda redshift,o_k,o_M,o_Lambda:1/((1+redshift)*np.sqrt(o_M*(1+redshift)**3+o_k*(1+redshift)**2+o_Lambda))
    tl = t_h * quad(E,0,z,args=(omega_k,omega_M,omega_Lambda))[0]
    return tl # in s

#convert m to kpc
def m_to_kpc(d):
    return d /(3.086*10**(16) *10**(3))

#convert m to Gpc
def m_to_gpc(d):
    return d /(3.086*10**(16) *10**(9))

#convert m to Gly
def m_to_gly(d):
    return d / 9.454254955488E+24

#convert Gpc to Gly
def gpc_to_gly(d):
    return d * 3.2638

#convert Gly to Gpc
def gly_to_gpc(d):
    return d / 3.2638

#convert s to Yr
def s_to_yr(t):
    return t/(60*60*24*365)

#convert s to GYr
def s_to_gyr(t):
    return t/(60*60*24*365*10**9)

def s_to_h(t,H0):
    return t*((H0/3.086)*10**(-19))

template = " |  {0:10} {1:29}|" # column widths: 10, 30
template2= " |  {0:10} {1:29}|" # column widths: 10, 30

print(" --------------World Model-------------------")
print(template.format(" "," "))
print(" | H : "+"%.1f" % H0+" km/s * Mpc       omega_m : "+"%.2f" % omega_M+" |")
print(" | omega_l : "+"%.2f" % omega_lambda+"            omega_k : "+"%.2f" % omega_k+" |")
print(" |                z : "+"%.1f" % z+"                   |")
print(template.format(" "," "))
print(" -------------Lookback-Time------------------")
print(template.format(" "," "))
print(template.format(str("%.3E" % t_l(omega_k,omega_M,omega_lambda,z,H0)),' seconds'))
print(template.format(str("%.3f" % s_to_h(t_l(omega_k,omega_M,omega_lambda,z,H0),H0)),' Hubble time'))
print(template.format(str("%.3f" % s_to_gyr(t_l(omega_k,omega_M,omega_lambda,z,H0))),' Giga years'))
print(template.format(str("%.3f" % s_to_gyr(1/(((H0/3.086)*10**(-19)))-t_l(omega_k,omega_M,omega_lambda,z,H0))),' Giga years after Big Bang'))
print(template.format(" "," "))

print(" -----------Comoving distance----------------")
print(template.format(" "," "))
print(template.format(str("%.3E" % D_M(omega_k,omega_M,omega_lambda,z,H0)),' meters'))
print(template.format(str("%.3f" % m_to_gpc(D_M(omega_k,omega_M,omega_lambda,z,H0))),' Giga parsecs'))
print(template.format(str("%.3f" % m_to_gly(D_M(omega_k,omega_M,omega_lambda,z,H0))),' Giga lightyears'))
print(template.format(" "," "))
print(template.format(str("%.3E" % m_to_kpc(np.pi/(180*60) * D_A(omega_k,omega_M,omega_lambda,z,H0))),' kpc/arcmin'))
print(template.format(str("%.3E" % m_to_kpc(np.pi/(180*3600) * D_A(omega_k,omega_M,omega_lambda,z,H0))),' kpc/arcsec'))
print(template.format(" "," "))
print(" --------------------------------------------")
