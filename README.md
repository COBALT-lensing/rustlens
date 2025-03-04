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