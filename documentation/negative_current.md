# Negative plasma current: assessment only

**No runtime polarity support is changed by this assessment.** The current reconstruction assumes
the existing positive-current, high-flux magnetic axis. Reversing current requires more than
negating the final result. Coordinate conventions remain COCOS 13; current polarity must not be
confused with changing coordinate convention.

## Linear parts

Green-function superposition, measurement replacement, the least-squares solve and the Poisson
operator are linear. Signed coil/passive/plasma currents can pass through these stages without
reversing R or Z coordinates. The fitted coefficients themselves may have either sign.
Regularisation penalties are quadratic and do not require a new sign.

The nonlinear axis, boundary and core decisions do assume a polarity:

| Stage | Present assumption | Required investigation/change |
|---|---|---|
| Magnetic axis | Choose the o-point with largest psi | Distinguish the plasma axis from other extrema; choose the appropriate maximum or minimum |
| X-point order | Descending psi | Order by an explicitly defined outward flux coordinate, not always raw descending psi |
| Limiter candidate | Largest visible limiter psi | Use the physically appropriate candidate in the same signed coordinate |
| Boundary selection | Axis psi exceeds boundary psi | Validate outward separation independent of raw flux sign |
| Flood fill | Core-side inequality assumes psi above the boundary | Use a sign-neutral inside/outside predicate, consistently with saddle handling |
| Normalisation | Ratio uses axis/boundary flux | The ratio can remain `(psi - axis_psi)/(boundary_psi - axis_psi)`; audit clamping, tolerances and rejection tests |
| Current centroid | Compute only for `ip > IP_MIN` | Test `abs(ip)` against a positive threshold, retain signed current in the integral ratio |
| Replay stop threshold | May use a signed current threshold | Define whether the threshold is magnitude or intentionally polarity-specific |
| Diagnostics | Field/contour signs and labels | Preserve signed fields and flux, positive geometry/penalties, document native-to-Wb conversion |

## Proposed design to evaluate

Introduce one explicitly defined flux orientation `s` and use `s * psi` only for topology
comparisons, not to overwrite the physical flux. Define `s` so the selected axis is a maximum in
the comparison coordinate. It must come from a documented current/flux convention, not an
unverified `sign(ip)` shortcut: COCOS, Green-function signs and source flux normalisation all
matter. Apply the same orientation to axis, saddle, limiter and core tests.

Zero current is not an orientation. Near the configured current floor, return an explicit
no-axis/invalid-equilibrium state rather than silently selecting a polarity. A reversal through
zero needs an agreed policy for retaining or resetting the previous Picard state.

The existing normalised flux and basis functions are promising because they are insensitive to
a simultaneous sign reversal of physical flux, axis flux and boundary flux. This does not by
itself prove that axis selection or flood filling is correct.

## Required acceptance tests before implementation is called complete

1. Replay a valid synthetic case, then reverse **all** contributing currents and magnetic signals
   with the same geometry. Verify signed psi, fields, plasma/passive currents and reconstructed
   signals reverse; axis location, normalised flux, mask, shape and chi-square stay consistent.
2. Cover limited and diverted equilibria, two X-points, competing o-points and saddles on grid
   vertices. Verify comparisons and tie-breaking, not just successful return codes.
3. Test near-zero current, exactly zero current and a time sequence crossing zero. Specify the
   expected state-reset/continuation policy before choosing a behaviour.
4. Verify C and offline Python conventions against an independent signed equilibrium and GSFit.
   Preserve the existing positive-current and v1 regression tolerances.
5. Audit external PCS assumptions and thresholds: fixing only the library is insufficient if
   the caller still rejects negative current.

The faster Poisson solver in #105 is independent work and remains deferred.
