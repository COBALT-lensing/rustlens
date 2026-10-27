# rustlens

A lightweight Python module for computing microlensing magnifications. rustlens is written in Rust.

## Installation

Install with pip:

```
pip install rustlens
```

## Usage

The plan is to eventually implement multiple lensing calculations. So far the only two available are based on [Witt & Mao (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract) and [Heyrovský (2003)](https://iopscience.iop.org/article/10.1086/376787/meta). Compute the magnification for uniform sources as follows:

```python
import rustlens

rustlens.witt_mao_magnification(l=[0.0, 0.2, 0.4, 0.6], rstar=2.0, re=0.5)
```

Where `l` is the distance of the lens from the centre of the source (as a list), normalised so that `l=1` is the edge of the disk, `rstar` is the source radius, and `re` is the Einstein ring radius (in any units as long as they're the same as each other).

To compute the brightness-integrated magnification for a non-uniform source, use:

```python
import rustlens

rustlens.integrated_witt_mao_magnification(l=[0.0, 0.2, 0.4, 0.6], rstar=2.0, re=0.5)
```

This assumes a linear limb-darkening law with coefficient 0.6.

You can also specify an arbitrary brightness/flux map:

```python
import rustlens

rustlens.integrated_flux_map_witt_mao_magnification(
    l=[0.0, 0.2, 0.4, 0.6],
    rstar=2.0,
    re=0.5,
    bl=[...],
    bf=[...],
)
```

Here `bl` contains strictly increasing radial coordinates in `[0, 1]`, and `bf` contains the corresponding nonnegative surface brightness. At least two points and positive total flux are required. The profile is linearly interpolated between these points and is zero outside their range. These are circularly symmetric radial profiles, not two-dimensional brightness maps.

It is possible to compute magnifications for a list of Einstein radii and stellar radii:

```python
import rustlens

rustlens.multi_witt_mao_magnification(l=[0, 0.5, 1], re=[0.1, 2, 5], rstar=[1, 2, 3])
# Returns: [[1.019803902718557, 1.0196620285602234, 1.00957643202704], [1.004987562112089, 1.0049780683135883, 1.0024469748497453], [1.002219758558194, 1.0022178584320351, 1.0010953955958362], [4.123105625617661, 3.881266049240038, 2.749075721239495], [2.23606797749979, 2.139190874964803, 1.6366197723675813], [1.6666666666666667, 1.6187416347356616, 1.3281528385688355], [10.04987562112089, 9.401134493155126, 6.450412687665968], [5.0990195135927845, 4.787617440381483, 3.3477740839861827], [3.4801021696368504, 3.2855490791489648, 2.360746191094734]]
```

## Heyrovský calculations

```python
rustlens.integrated_heyrovsky_magnification(l=[0.0, 0.5, 1.0], rstar=2.0, re=0.5)
rustlens.integrated_flux_map_heyrovsky_magnification(
    l=[0.0, 0.5, 1.0], rstar=2.0, re=0.5, bl=[0.0, 1.0], bf=[1.0, 1.0]
)
```

These return total magnification, normalized by unlensed flux, with the same
brightness conventions as the Witt–Mao functions.

`heyrovsky_magnification(l, r, re, rstar)` instead returns the angular mean
magnification of an infinitesimally thin ring at source-normalized radius `r`.
It is the paper's angular integral divided by `2*pi`, not a disk magnification.
It returns infinity at `l == r`; integrating a bounded surface-brightness profile
over the disk gives a finite magnification.

For `discrete_flux_map_heyrovsky_magnification`, `bf` must contain unlensed
**ring flux weights**, including area and quadrature weights, rather than surface
brightness samples. The result is `sum(bf * ring_magnification) / sum(bf)`.
Previous versions omitted this normalization. Use the integrated flux-map
function for sampled surface brightness.

Radii must be positive and finite, in the same projected plane (or both angular);
separations must be finite and nonnegative. Invalid inputs raise `ValueError`;
failed numerical evaluations raise `RuntimeError`.

See [the mathematical audit](MATHEMATICS.md) for equation references, corrections,
numerical treatment, and verification limits. Run the regression checks with
`cargo test --lib`.
