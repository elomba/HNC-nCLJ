"""Classical equal-partial-volume mixture 2PT, with signed R2PT bookkeeping.

References: Sun et al., JCP 147,194505 (2017), Eqs 2,3,8-12,A3-A7;
Lai et al., PCCP 14,15206 (2012), equal-partial-volume one-fluid approximation.
No parameters in this module are fitted to thermodynamic-integration entropy.
All particles have unit mass, kB=1. Frequency nu is cycles/time, not radians/time.
"""
import numpy as np
from scipy.fft import dct, next_fast_len
from scipy.optimize import brentq


def trap(y,x):
    return float(np.trapezoid(y,x))


def delta_from_D(D,T,rho_eff,mass=1.):
    # Sun Eq A3; crucially rho^(+1/3), not rho^(-1/3).
    return (8/3)*(36/np.pi**2)**(1/3)*D*np.sqrt(np.pi*mass/T)*rho_eff**(1/3)


def fluidicity(delta,exponent=1.):
    # Solve Sun A4/A7 in log(gamma) form, avoiding large polynomial powers.
    assert delta>0 and exponent>0
    def equation(log_gamma):
        gamma=np.exp(log_gamma)
        logf=(log_gamma+1.5*np.log(delta))/(1+1.5*exponent)
        return (2/3)*logf-(np.log(delta/2)+(2/3)*log_gamma+np.log(2-gamma)-3*np.log1p(-gamma))
    logg=brentq(equation,-200.,-1e-12,xtol=1e-13)
    gamma=np.exp(logg)
    f=np.exp((logg+1.5*np.log(delta))/(1+1.5*exponent))
    assert 0<f<1 and 0<gamma<1
    return float(f),float(gamma)


def log_moment(nu,G):
    """Integral G(nu)*ln(nu) using piecewise-linear G, including nu=0 cell.

    No arbitrary low-frequency epsilon; integrate logarithmic singularity exactly.
    """
    x=np.asarray(nu);g=np.asarray(G)
    b=np.diff(g)/np.diff(x);a=g[:-1]-b*x[:-1]
    xl=np.zeros_like(x);ok=x>0;xl[ok]=x[ok]*np.log(x[ok])
    F1=xl-x;F2=.5*x*xl-.25*x*x
    return float(np.sum(a*np.diff(F1)+b*np.diff(F2)))


def spectrum(c_biased,dt,taper_tau,T,padding_factor=1):
    """Nonnegative Gaussian-smoothed periodogram via full biased autocorrelation.

    Biased C(lag)=sum(v(t+lag).v(t))/(3*Nsample*Nframes) is an exact
    finite-record autocorrelation, so its Fourier transform is nonnegative.
    An optional Gaussian multiplier preserves positive definiteness. We retain
    every lag and the complete Nyquist band, rather than truncating/renormalizing
    a narrow frequency band. C0 normalization enforces equipartition; its raw
    target-T mode count and kinetic discrepancy are explicitly returned.
    """
    c=np.asarray(c_biased,float);C0=float(c[0]);assert C0>0
    time=np.arange(len(c))*dt
    taper=np.ones_like(time) if taper_tau is None else np.exp(-(time/taper_tau)**2)
    size=next_fast_len(padding_factor*len(c))+1
    # size must allow the autocorrelation to become zero outside the record.
    if size<=len(c):size=next_fast_len(2*len(c)-1)+1
    a=np.zeros(size);a[:len(c)]=c*taper
    G=6*dt/C0*dct(a,type=1)
    nu=np.arange(size)/(2*(size-1)*dt)
    modes=trap(G,nu)
    assert abs(modes-3)<1e-10,(modes,dt,len(c))
    assert G.min()>-1e-8,float(G.min())
    D_raw=dt*(np.sum(c*taper)-.5*C0)  # exact same trapezoid as zero-frequency transform
    D_equipped=D_raw*T/C0
    assert abs(G[0]-12*D_equipped/T)<1e-10
    return nu,G,dict(C0=C0,raw_targetT_modes=3*C0/T,normalized_modes=modes,
                    D_spectral_raw=float(D_raw),D_spectral_equipartition=float(D_equipped),
                    s0=float(G[0]),min_spectrum=float(G.min()),
                    nyquist_frequency=float(nu[-1]),upper_10pct_mode_fraction=trap(G[nu>=.9*nu[-1]],nu[nu>=.9*nu[-1]])/3,
                    taper_tau=taper_tau)


def entropy_from_moments(modes,logmoment,S0,T,rho,rho_eff,exponent=1.,h=1.,hs_lnZ=False):
    """Absolute classical 2PT minus ideal mixture reference (mixing cancels).

    Gas entropy is f*s_HS(T,rho_eff*f), not the unweighted HS entropy.
    Solid entropy uses signed G-Ggas and the analytic Lorentzian tail to infinity.
    The arbitrary h cancels because the total mode weight is 3.
    hs_lnZ is an explicitly nonprimary historical same-P vs same-rho diagnostic.
    """
    assert S0>0
    D=S0*T/12
    delta=delta_from_D(D,T,rho_eff)
    f,gamma=fluidicity(delta,exponent)
    a=np.pi*S0/(6*f)
    sid=2.5+1.5*np.log(2*np.pi*T/h**2)-np.log(rho)
    shs=gamma*(3*gamma-4)/(1-gamma)**2
    Z=(1+gamma+gamma**2-gamma**3)/(1-gamma)**3
    sg=f*(2.5+1.5*np.log(2*np.pi*T/h**2)-np.log(rho_eff*f)+shs+(np.log(Z) if hs_lnZ else 0.))
    I_total=modes*(1+np.log(T/h))-logmoment
    I_gas=3*f*(1+np.log(a*T/h))
    ss=I_total-I_gas
    return dict(s_ex=float(sg+ss-sid),s_abs=float(sg+ss),s_id_reference=float(sid),
                s_gas_absolute=float(sg),s_solid_signed_absolute=float(ss),f=f,gamma=gamma,
                delta=float(delta),rho_eff=float(rho_eff),hs_excess_per_gas_particle=float(shs),
                hs_lnZ=hs_lnZ,exponent=exponent)


def evaluate(nu,G,T=2.,rho=.7,rho_eff=None,exponent=1.,h=1.,diagnostics=True,hs_lnZ=False):
    if rho_eff is None:rho_eff=rho
    modes=trap(G,nu);lm=log_moment(nu,G)
    out=entropy_from_moments(modes,lm,float(G[0]),T,rho,rho_eff,exponent,h,hs_lnZ)
    if diagnostics:
        f=out['f'];a=np.pi*G[0]/(6*f)
        gas=G[0]/(1+(a*nu)**2);residual=G-gas
        tail_modes=3*f-trap(gas,nu)
        negative=trap(np.maximum(-residual,0),nu)+tail_modes
        clipped=np.maximum(residual,0)
        clipped_entropy=trap(clipped,nu)*(1+np.log(T/h))-log_moment(nu,clipped)
        out.update(total_modes=modes,gas_modes_exact=3*f,solid_modes_exact=modes-3*f,
                   gas_tail_modes_above_nyquist=tail_modes,negative_residual_modes=negative,
                   negative_residual_fraction=negative/3,
                   clipped_total_mode_excess=3*f+trap(clipped,nu)-3,
                   clipped_s_ex=float(out['s_gas_absolute']+clipped_entropy-out['s_id_reference']),
                   spectrum_reconstruction_error=float(np.max(abs(G-(gas+residual)))))
    return out
