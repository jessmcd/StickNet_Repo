# Various thermodynamics functions used in plotting Skew-Ts

import numpy as np 

C_to_K = 273.15 
c_p_dry = 1005.7 # in J/kgK +- 2.5
c_V_dry = c_p_dry - 287.04 # Based on the given value for R-dry in J/kgK
eps = 0.6220
k_dry = 0.2854 # +- 0.0007

def sat_vapor_pressure(T):
	""" Saturation vapor pressure (mb) based on temperature in degrees Celsius from Bolton (1980) eqn 10 """
	es_t = 6.112 * np.exp(17.67*T/(T+243.5))
	return es_t

def sat_vapor_temperature(e_s):
	""" Calculated the temperature (in C) for a specific saturation vapor pressure (mb) from Bolton (1980) eqn 11"""
	T = (243.5 * np.log(e_s) - 440.8)/(19.48 - np.log(e_s))
	return T

def sat_mixing_ratio(p, T):
	""" Calculates the saturation mixing ratio from p (in mb) and T (in C)
	 """
	w_s = eps * sat_vapor_pressure(T) / (p - sat_vapor_pressure(T) )	
	return w_s

def mixing_ratio_line(p, w_s):
	""" Given pressure (mb) and saturation mixing ratio (unitless) calculate temperature (in C)  calculated from w=eps*e_s/(p-e_s)"""	
	e_s = w_s * p / (eps + w_s)
	T = sat_vapor_temperature(e_s)
	return T

def RH(T, p, w):
	""" Given the temperature (in C), the pressure (mb), and the mixing ratio (unitless) calculates the relative humidity (%)
	given RH = 100 * e / e_s"""
	DT = mixing_ratio_line(p, w)
	RH_p = 100. * sat_vapor_pressure(DT) / sat_vapor_pressure(T)
	return RH_p

def T_LCL(T, RH):
	""" Find the LCL temperature (K) given T in (K) and RH (%) from Bolton (1980) eqn 22"""
	TL = 1. / ((1./(T - 55.)) - (np.log(RH/100.)/2840.)) + 55.
	return TL

def theta_dry(theta, p, p_O = 1000.0):
	""" Find the dry potential temperature in K given the potential temperature and pressure(s) (in mb) where e=0 using Bolton (1980) eqn 23 (theta_dry = T_K * (p_0/(p-e))^k_dry) """
	T_K = theta*((p_O/p))**(-k_dry)
	return T_K

def pseudoeq_potential_T(T,p,w,p_O=1000.0):
	""" Calculates the pseudoadiabatic equivalent potential temperature (K) given T (C), p (mb), w (kg/kg), by Bolton (1980) eqn 43"""
	theta_ep = (T+C_to_K) * (p_O/p)**(0.2854*(1-0.28*w)) * np.exp(((3.376/T_LCL(T+C_to_K, RH(T, p, w))) - 0.00254) * w * 10**(3) *(1+0.81*w))
	return theta_ep

def theta_ep_field(T, p, p_O=1000.0):
	""" Given T (in C) and p (in mb) returns Theta_ep in K for a theoretical saturated parcel """
	w = sat_mixing_ratio(p, T)
	theta_ep_field = pseudoeq_potential_T(T, p, w)
	return theta_ep_field
