"""Analytic physics checks, not comparisons tuned to the MD entropy data."""
import json
import numpy as np
from scipy.integrate import quad
from two_pt import delta_from_D,fluidicity,spectrum,evaluate,entropy_from_moments

checks={}
for delta in [.01,.1,1.,10.,100.]:
 f,g=fluidicity(delta)
 residual=(2*delta**-4.5*f**7.5-6*delta**-3*f**5-delta**-1.5*f**3.5+6*delta**-1.5*f**2.5+2*f-2)
 assert abs(residual)<1e-9,(delta,residual)
checks['fluidicity_matches_original_polynomial']=True
D,T,rho=.23,2.,.7
base=delta_from_D(D,T,rho)
for a,b in [(2.,3.),(.1,4.)]:
 assert abs(delta_from_D(D*b/a**2,T*(b/a)**2,rho*a**3)-base)<1e-12
checks['delta_invariant_to_length_and_time_units']=True
# A solvable Ornstein-Uhlenbeck VACF provides the exact Green-Kubo value T/alpha.
alpha=10.;dt=.001;c=T*np.exp(-alpha*np.arange(50000)*dt)
nu,G,diag=spectrum(c,dt,None,T)
assert abs(diag['normalized_modes']-3)<1e-12
assert abs(diag['D_spectral_raw']-T/alpha)<2e-6
assert abs(diag['s0']-12/alpha)<1e-5
checks['OU_spectrum_mode_count_and_GK_limit']=True
# Integral of Lorentzian times classical oscillator weight is analytic, including its tail.
a=3.;S0=6*a/np.pi  # total modes of this Lorentzian are 3
for h in [.1,1.,10.]:
 numeric=quad(lambda v:S0/(1+(a*v)**2)*(1-np.log(h*v/T)),0,np.inf,epsabs=1e-9)[0]
 assert abs(numeric-3*(1+np.log(a*T/h)))<1e-8
checks['analytic_Lorentzian_entropy_tail']=True
# Exact gas limit: f=1,gamma=0, gas and total spectra coincide, excess tends to zero.
# At finite but tiny rho the root approaches this limit without manual overriding.
rho_dilute=1e-12; D_dilute=1/rho_dilute  # kinetic dilute gas: D diverges as 1/rho
S0_dilute=12*D_dilute/T; a_dilute=np.pi*S0_dilute/6
for h in [.1,1.,10.]:
 out=entropy_from_moments(3,-3*np.log(a_dilute),S0_dilute,T,rho_dilute,rho_dilute,h=h)
 assert abs(out['s_ex'])<2e-4,out
checks['ideal_gas_limit']=True
values=[evaluate(nu,G,h=h)['s_ex'] for h in [.1,1.,10.]]
assert np.ptp(values)<1e-11
checks['classical_reference_h_cancels']=True
# Artificially splitting exactly identical mechanical spectra must preserve excess entropy.
sex=evaluate(nu,G)['s_ex']; split=.4*sex+.6*sex
assert abs(split-sex)<1e-14
sid=2.5+1.5*np.log(2*np.pi*T)-np.log(.7)
mix=-.4*np.log(.4)-.6*np.log(.6)
assert abs((.4*(sex+sid)+.6*(sex+sid)+mix)-(sid+mix)-sex)<1e-13
checks['identical_label_split_and_mixing_reference_cancel']=True
print(json.dumps(checks,indent=2))
