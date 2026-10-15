# Solution explanation — seasonal dengue transmission potential

Task package: `tasks/biology/mathematical-epidemiology/seasonal-dengue-transmission-potential`

## The research question and why it is hard

Mechanistic temperature-dependent models of *Aedes*-borne disease are routinely used by public-health
modelling groups to decide whether a city is at risk of autochthonous dengue transmission after an
importation. The workflow has three layers that must be done in order. First, laboratory thermal
performance curves for mosquito and virus traits are estimated. Second, they are combined with the local
climate and a field estimate of mosquito abundance. Third, a threshold quantity, the basic reproduction
number, is computed for the seasonal environment. Each layer feeds the next. An error in temperature
units or in the abundance calibration propagates multiplicatively into R0. Several plausible shortcuts
for the last layer also give different answers on this dataset: R0 at the mean temperature gives about
1.07, the time-averaged instantaneous R0 gives 1.78, and a quasi-static mosquito density gives 1.72,
against a correct value of about 1.53. Two of those shortcuts sit far from the threshold and one sits
right on it, so the method changes the public-health conclusion.

## Interpreting and harmonising the data

An expert starts by reading `lab_registry.csv` alongside the trait table, because the three laboratories
recorded their experiments differently. FL-VBL reported temperatures in °F and probabilities as
percentages. It also summarised adult survival as a constant daily mortality rate, and under a constant
hazard the mean lifespan is lf = 1/μ. GZ-INS reported the extrinsic incubation period as the time to 50 %
dissemination; following the usual convention in thermal-trait syntheses, the rate in the model is
PDR = 1/EIP. KH-ENT recorded controller setpoints rather than realised chamber temperatures, and an audit
documents a uniform +1.5 °C offset. Ignoring that offset shifts every KH-ENT trait by 1.5 °C. The worst
case is egg production near its upper thermal limit, where the fitted Tm moves enough to change R0 at
the mean temperature by about 16 %. Replicates flagged by QA for chamber failure contain values
inflated or deflated several-fold and must be excluded. Least squares is not robust to them, and keeping
them distorts the development-rate curve by about 30 %. Zeros at extreme temperatures are real
observations (no development, no transmission) and must stay in, because they identify T0 and Tm.

## Fitting the thermal performance curves

Each trait follows the Brière or quadratic form prescribed in the specification. The reference solution
minimises the residual sum of squares over all retained replicates. Because q enters linearly once
(T0, Tm) are fixed, the solution first profiles q out in closed form over a grid of (T0, Tm). It then
refines all three parameters with Nelder–Mead under the physical constraints T0 < min{T : y > 0} and
Tm > max{T : y > 0}. This avoids the local minima the piecewise-defined curves create for gradient
methods started far from the optimum. Other reasonable estimators, such as fits to temperature means,
binomial likelihood for the probabilities, or log-scale least squares for the multiplicative rate noise,
give curves within a few per cent in the data-supported range, and the verifier tolerates them.

## Climate

The station series needs three corrections before it represents the urban area: removal of the
missing-data code and blanks, homogenisation of the pre-2018 hilltop record (+1.1 °C, as stated in the
station history), and removal of the faulty-probe episode. The maintenance log only says the sensor was
replaced on 21 September 2016 and that the fault onset is unknown, so the expert has to find the fault in
the data: a 49-day run of the identical value 18.0 °C in late summer. Leaving it in lowers the
climatology by about 0.15 °C overall, and by much more in August–September when transmission peaks. The
periodic forcing T(t) is then the mean annual cycle, estimated by least squares on a few annual harmonics
after rescaling leap years to 365 days. A smoothed day-of-year climatology is an equally valid choice and
gives the same answers within tolerance.

## Disease-free seasonal vector state and calibration of δ

The immature–adult system is cooperative and its nonlinearity is strictly sublinear. For such a
periodic system, if the zero solution is unstable there is a unique, globally attractive positive
periodic solution (Zhao, 2017). It is computed by integrating over successive years until the state at
the start of the year stops changing. One structural observation makes the unknown δ easy to handle.
Substituting L = ℓ/δ and M = m/δ removes δ from the equations, so M*(t) = m*(t)/δ exactly. The abundance
scale k = 1/(δ N_h) therefore enters linearly, and it can be estimated by maximum likelihood from the
survey. Treating weekly counts as Poisson with mean 0.55 × residents × k × (weekly mean of m*), the
estimate is k̂ = ΣC / Σ(0.55 R m̄). The 55 % aspirator efficiency matters. Without it, density is underestimated by a factor
of 0.55 and R0 by a factor of √0.55. Because the population lags temperature by about 60 days (peak around day 256 against a
temperature maximum near day 198), calibrating with the quasi-equilibrium m_eq(T(t)) instead of m*(t)
gives the wrong k and the wrong timing. The immature equation becomes stiff in summer through the −δL²
term. An explicit fixed-step RK4 scheme with an eighth-of-a-day step is stable but inaccurate there,
biasing the summer peak by about 15 % and k by about 10 %. An implicit or adaptive solver is needed.

## Reproduction numbers

At constant temperature the next-generation matrix has only cross terms. Its spectral radius is
R0(T) = sqrt{ [a b PDR / (μ(PDR + μ))] · [a c k m_eq(T) / γ] } with μ = 1/lf. The first factor is the
number of humans infected by one newly infected mosquito, with survival through the extrinsic incubation
period. The second is the number of mosquitoes infected by one human over the viraemic period. The
equilibrium density is m_eq = ½ MDR² (½ EFD/μ − 1/pEA)/μ, which is positive only when the mosquito net
reproduction number ½ EFD·pEA/μ exceeds 1. For the seasonal climate, R0 is defined as the spectral
radius of the next-generation operator on 365-periodic functions (Bacaër & Guernaoui, 2006). By Wang and
Zhao (2008, Theorem 2.2), R0 is the unique λ > 0 for which the monodromy matrix of
x′ = (F(t)/λ − V(t)) x has spectral radius 1. The theorem applies here because the transmission-free
linear system x′ = −V(t)x is exponentially stable: all its diagonal rates are positive. The solution
computes the monodromy matrix with monthly renormalisation to avoid overflow and solves for λ with
Brent's method. A day-step discretisation of the operator converges to the same value. The
time-averaged instantaneous R0 and R0 at the mean temperature are not threshold quantities for the
periodic system. The first ignores the seasonal lag between temperature, mosquito density and
extrinsic incubation; the second ignores the strong convexity of R0(T) between 17 and 29 °C.

## How correctness is verified

The verifier recomputes everything from the raw files with its own implementation, using different
numerics: multi-start bounded least squares, an implicit Radau integrator for the vector state, and
fixed-step RK4 with bisection for the Floquet problem. It then checks four things. (i) Each fitted
thermal curve against the independent fit over 18–34 °C. (ii) Internal consistency of the
constant-temperature table with the agent's own fitted traits: the density column must be a constant
multiple of the specified equilibrium, and R0 must follow the formula. (iii) Every summary quantity
against the reference, with tolerances set between the spread of valid methods (under 1 % for seasonal
R0) and the effect of each documented error (4–16 %). (iv) Seasonal R0 recomputed from the agent's own
traits and density scale, so a correct periodic computation built on slightly different fits is
accepted. Eleven deliberately broken pipelines, one per pitfall described above, all fail. The oracle and
two alternative pipelines that use different valid estimators all pass.
