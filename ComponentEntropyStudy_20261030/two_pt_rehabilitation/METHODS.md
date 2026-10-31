# Classical 2PT formulation used in this test

The implementation is `two_pt.py`. It follows the signed-residual construction of [Sun et al. (2017)](https://doi.org/10.1063/1.5001798), with an explicitly approximate mixture-volume prescription from [Lai et al. (2012)](https://doi.org/10.1039/C2CP42011B). Their excess Gibbs energy of mixing is a different reference from the ideal-gas excess entropy used in the present Rosenfeld test; the volume prescription is borrowed, not their numerical validation.

All particles have mass one and kB=1. Frequency ν is cycles per LJ time. For each distinct mechanical species, C(t)=mean[v(t)·v(0)]/3 after subtracting the same whole-system COM velocity. The finite-record autocorrelation uses a fixed denominator equal to the number of frames, giving a nonnegative periodogram. No particle-wise temporal mean is subtracted. The Gaussian-tapered spectrum is

\[
G(\nu)=\frac{12}{C(0)}\int_0^\infty C_{\rm biased}(t)e^{-(t/\tau_w)^2}\cos(2\pi\nu t)\,dt.
\]

The full discrete Nyquist band is retained, with exactly three modes after equipartition normalization. The unnormalized target-temperature mode count is saved separately. The zero-frequency value gives the corresponding smoothed, kinetic-normalized diffusion estimate D=G(0)T/12; it is compared with the independently fitted MSD coefficient, not assumed equal to it.

The equal-partial-volume primary approximation is Vbar_i=1/rho. Thus rho_eff=1/Vbar_i=rho; physical composition weights remain xi. Different actual partial volumes are not ruled out by equal nominal LJ diameters. Volume-conserving variations are reported as a stress test. Species allocations here are not independently evaluated partial molar derivatives.

The dimensionless diffusivity and partition equations are

\[
\Delta=\frac83\left(\frac{36}{\pi^2}\right)^{1/3}
 D\sqrt{\frac{\pi}{T}}\,\rho_{\rm eff}^{1/3},\qquad
\gamma=\Delta^{-3/2}f^{1+3\delta/2},
\]
\[
f^{2/3}=\frac\Delta2\gamma^{2/3}\frac{2-\gamma}{(1-\gamma)^3}.
\]

They are solved in log(γ) for physical 0<f,γ<1. The primary exponent δ=1 recovers the standard fluidicity polynomial. The published LJ sensitivity δ=1.5 is not refitted to these data.

For a=πG(0)/(6f), the gas spectrum is Gg(ν)=G(0)/[1+(aν)²], and its integral to infinity is 3f. The classical oscillator weight is Ws(ν)=1−ln(hν/T). With s_id(rho)=5/2+(3/2)ln(2πT/h²)−ln(rho),

\[
s_g=f\left[\frac52+\frac32\ln\frac{2\pi T}{h^2}
-\ln(f\rho_{\rm eff})+\frac{\gamma(3\gamma-4)}{(1-\gamma)^2}\right],
\]
\[
s_s=\int_0^\infty G(\nu)W_s(\nu)d\nu
-3f\left[1+\ln\frac{aT}{h}\right],\qquad
s_{{\rm ex},i}=s_g+s_s-s_{\rm id}(\rho).
\]

The second term in ss is the analytical gas-spectrum oscillator integral, including the infinite Lorentzian tail. The total-spectrum logarithmic moment is integrated using piecewise-linear spectra, including the ν=0 cell analytically. This avoids a low-frequency epsilon. No residual is clipped: G−Gg can be negative as bookkeeping in the revised model, not as a literal population of physical oscillators.

The mixture entropy estimate is Σxi s_ex,i. Adding ideal mixing −Σxi ln xi to the approximate absolute mixture entropy and subtracting it in the ideal-gas reference cancels it exactly. Because the total spectral weight is three, h also cancels. Both cancellations are tested. Mechanically identical labels are merged with matched retained IDs before the nonlinear spectral transformation; separately processed labels are retained as a finite-sample diagnostic.

The primary smoothing width is 10 tau. Other smoothing widths, record halves, frequency-grid refinement, sampling decimation, partition exponent, and partial-volume perturbations are sensitivities. None is selected by agreement with TI. The 0.05 kB absolute-bias and 0.02 kB increment-bias tolerances in `protocol.json` are practical thresholds for this study, not universal standards. A point estimate within tolerance and an interval containing zero are consistency checks, not an equivalence proof at that tolerance.
