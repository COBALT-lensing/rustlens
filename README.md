# rustlens

A lightweight Python module for computing microlensing magnifications. rustlens is written in Rust and has no external dependencies.

## Installation

Install with pip:

```
pip install rustlens
```

## Usage

The plan is to eventually implement multiple lensing calculations. So far the only one available is based on [Witt & Mao (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract). Compute the magnification for uniform sources as follows:

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

This will assume the source is limb darkened using a linear limb darkening law. 

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

Where `bl` is a list of coordinates and `bf` is the corresponding brightness at each point. These values will be interpolated to integrate the brightness profile.

It is possible to compute magnifications for a list of Einstein radii and stellar radii:

```python
import rustlens

rustlens.multi_witt_mao_magnification(l=[0, 0.5, 1], re=[0.1, 2, 5], rstar=[1, 2, 3])
# Returns: [[1.019803902718557, 1.0196620285602234, 1.00957643202704], [1.004987562112089, 1.0049780683135883, 1.0024469748497453], [1.002219758558194, 1.0022178584320351, 1.0010953955958362], [4.123105625617661, 3.881266049240038, 2.749075721239495], [2.23606797749979, 2.139190874964803, 1.6366197723675813], [1.6666666666666667, 1.6187416347356616, 1.3281528385688355], [10.04987562112089, 9.401134493155126, 6.450412687665968], [5.0990195135927845, 4.787617440381483, 3.3477740839861827], [3.4801021696368504, 3.2855490791489648, 2.360746191094734]]
```