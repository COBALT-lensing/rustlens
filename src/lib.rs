use core::f64;
use interp::interp;
use interp::InterpMode;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use quadrature::integrate;

const LD_COEFF: f64 = 0.6;

fn ld_linear(r: f64) -> f64 {
    1.0 - LD_COEFF * (1.0 - (1.0 - r.powi(2)).sqrt())
}

// Public coordinates use the source radius; elliptic disk coordinates use re.
fn validate(l: &[f64], re: f64, rstar: f64) -> PyResult<()> {
    if !re.is_finite()
        || re <= 0.0
        || !rstar.is_finite()
        || rstar <= 0.0
        || !(rstar / re).is_finite()
        || !(re / rstar).is_finite()
        || re / rstar == 0.0
        || l.iter().any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(PyValueError::new_err(
            "radii must be positive and finite; separations must be finite and nonnegative",
        ));
    }
    Ok(())
}

// quadrature's fixed refinement budget can be exhausted, and it silently drops
// non-finite samples. Detect both conditions instead of accepting a false result.
fn quad(f: &dyn Fn(f64) -> f64, a: f64, b: f64, tol: f64, depth: u32) -> f64 {
    if a == b {
        return 0.0;
    }
    let bad = std::cell::Cell::new(false);
    let q = integrate(
        |x| {
            let y = f(x);
            if !y.is_finite() {
                bad.set(true);
            }
            y
        },
        a,
        b,
        tol,
    );
    if bad.get() {
        return f64::NAN;
    }
    if q.error_estimate <= tol.max(2e-12 * q.integral.abs()) {
        return q.integral;
    }
    if depth == 0 {
        return f64::NAN;
    }
    let mid = a + (b - a) / 2.0;
    quad(f, a, mid, tol / 2.0, depth - 1) + quad(f, mid, b, tol / 2.0, depth - 1)
}

fn checked(value: f64) -> PyResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PyRuntimeError::new_err(
            "magnification did not converge to a finite value",
        ))
    }
}

fn profile(bl: &[f64], bf: &[f64]) -> PyResult<f64> {
    if bl.len() < 2
        || bl.len() != bf.len()
        || bl.iter().any(|r| !r.is_finite() || *r < 0.0 || *r > 1.0)
        || bl.windows(2).any(|w| w[0] >= w[1])
        || bf.iter().any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(PyValueError::new_err("profile needs matching arrays of at least two strictly increasing radii in [0, 1] and finite nonnegative brightnesses"));
    }
    // Exact integral of r B(r) for each linear segment; outside the knots B=0.
    let flux: f64 = bl
        .windows(2)
        .zip(bf.windows(2))
        .map(|(r, b)| {
            (r[1] - r[0]) / 6.0 * ((2.0 * r[0] + r[1]) * b[0] + (r[0] + 2.0 * r[1]) * b[1])
        })
        .sum();
    if !flux.is_finite() || flux <= 0.0 {
        return Err(PyValueError::new_err(
            "profile must have positive finite unlensed flux",
        ));
    }
    Ok(flux)
}

#[pyfunction]
pub fn integrated_heyrovsky_magnification(l: Vec<f64>, re: f64, rstar: f64) -> PyResult<Vec<f64>> {
    validate(&l, re, rstar)?;
    l.iter()
        .map(|&x| {
            integrate_rings(
                x,
                re / rstar,
                &ld_linear,
                &[0.0, 1.0],
                (1.0 - LD_COEFF / 3.0) / 2.0,
            )
        })
        .collect()
}

#[pyfunction]
pub fn discrete_flux_map_heyrovsky_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    bl: Vec<f64>,
    bf: Vec<f64>,
) -> PyResult<Vec<f64>> {
    validate(&l, re, rstar)?;
    if bl.len() != bf.len()
        || bl.is_empty()
        || bl.iter().any(|r| !r.is_finite() || *r < 0.0 || *r > 1.0)
        || bf.iter().any(|f| !f.is_finite() || *f < 0.0)
    {
        return Err(PyValueError::new_err(
            "discrete rings need matching radii in [0, 1] and nonnegative finite flux weights",
        ));
    }
    let total: f64 = bf.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(PyValueError::new_err(
            "flux weights must have a positive finite sum",
        ));
    }
    // bf is unlensed flux per infinitesimal ring (including area/quadrature weights).
    // An exact ring-lens coincidence really does have infinite magnification.
    l.iter()
        .map(|&x| {
            let mut result = 0.0;
            for (&r, &f) in bl.iter().zip(&bf) {
                if f == 0.0 {
                    continue;
                }
                if x == r {
                    return Ok(f64::INFINITY);
                }
                result += checked(ring_magnification(x, r, re / rstar))? * (f / total);
            }
            checked(result)
        })
        .collect()
}

#[pyfunction]
pub fn integrated_flux_map_heyrovsky_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    bl: Vec<f64>,
    bf: Vec<f64>,
) -> PyResult<Vec<f64>> {
    validate(&l, re, rstar)?;
    let flux = profile(&bl, &bf)?;
    let b = |r| interp(&bl, &bf, r, &InterpMode::Constant(0.0));
    l.iter()
        .map(|&x| integrate_rings(x, re / rstar, &b, &bl, flux))
        .collect()
}

fn integrate_rings(
    l: f64,
    epsilon: f64,
    b: &dyn Fn(f64) -> f64,
    knots: &[f64],
    flux: f64,
) -> PyResult<f64> {
    let mut points = knots.to_vec();
    if l > knots[0] && l < knots[knots.len() - 1] {
        points.push(l);
    }
    points.sort_by(f64::total_cmp);
    points.dedup();
    let mut sum = 0.0;
    for w in points.windows(2) {
        let (a, z) = (w[0], w[1]);
        let weighted = |r: f64| {
            if r == 0.0 {
                return if l == 0.0 {
                    epsilon * b(0.0) / flux
                } else {
                    0.0
                };
            }
            r * b(r) / flux * ring_magnification(l, r, epsilon)
        };
        // Move the logarithmic singularity to t=0 and cancel it with dr/dt.
        let value = if a == l || z == l {
            quad(
                &|t| {
                    let r = if a == l {
                        a + (z - a) * t * t
                    } else {
                        z - (z - a) * t * t
                    };
                    if t == 0.0 || r == l {
                        return 0.0;
                    }
                    2.0 * (z - a) * t * weighted(r)
                },
                0.0,
                1.0,
                1e-10 / points.len() as f64,
                12,
            )
        } else {
            quad(&weighted, a, z, 1e-10 / points.len() as f64, 12)
        };
        sum += value;
    }
    checked(sum)
}

#[pyfunction]
pub fn heyrovsky_magnification(l: f64, r: f64, re: f64, rstar: f64) -> PyResult<f64> {
    validate(&[l, r], re, rstar)?;
    let result = ring_magnification(l, r, re / rstar);
    if result.is_infinite() && l == r {
        return Ok(result);
    }
    checked(result)
}

fn ring_magnification(l: f64, r: f64, epsilon: f64) -> f64 {
    if l == r {
        return f64::INFINITY;
    }
    if l == 0.0 || r == 0.0 {
        return point_source_magnification(l.max(r) / epsilon);
    }
    // Heyrovsky (2003), eq. (4), divided by 2 pi: the angular MEAN.
    // Evaluate the same K/Pi combination in Bulirsch form so that 1-n and
    // 1-m do not lose all precision near r=l. No asymptotic cutoff is needed.
    let scale = l.max(r).max(epsilon);
    let d = (l - r).abs() / scale;
    let (l, r, e) = (l / scale, r / scale, epsilon / scale);
    let s = l + r;
    let h = d.hypot(2.0 * e);
    let p = (d / s).powi(2);
    let kc = (d / s) * (s.hypot(2.0 * e) / h);
    let a = 2.0 * e * e + d * d;
    let b = 2.0 * e * e * p + d * d;
    match ellip::cel(kc, p, a, b) {
        Ok(v) => 2.0 * v / (f64::consts::PI * s * h),
        Err(_) => f64::NAN,
    }
}

#[pyfunction]
pub fn integrated_witt_mao_magnification(l: Vec<f64>, re: f64, rstar: f64) -> PyResult<Vec<f64>> {
    validate(&l, re, rstar)?;
    let rho = rstar / re;
    l.iter()
        .map(|&x| {
            let c = |r: f64| {
                if r == 0.0 {
                    0.0
                } else {
                    r * r * witt_mao_magnification_scaled(x * rho, r * rho)
                }
            };
            // Integration by parts of integral B(r) d[r^2 A_disk(r)].
            // mu=sqrt(1-r^2) removes the derivative's limb singularity.
            let f = |mu: f64| c((1.0 - mu * mu).max(0.0).sqrt());
            let cut = if x < 1.0 { (1.0 - x * x).sqrt() } else { 0.0 };
            let integral = quad(&f, 0.0, cut, 1e-10, 12) + quad(&f, cut, 1.0, 1e-10, 12);
            checked(((1.0 - LD_COEFF) * c(1.0) + LD_COEFF * integral) / (1.0 - LD_COEFF / 3.0))
        })
        .collect()
}

#[pyfunction]
pub fn integrated_flux_map_witt_mao_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    bl: Vec<f64>,
    bf: Vec<f64>,
) -> PyResult<Vec<f64>> {
    validate(&l, re, rstar)?;
    let flux = profile(&bl, &bf)?;
    let rho = rstar / re;
    l.iter()
        .map(|&x| {
            let c = |r: f64| {
                if r == 0.0 {
                    0.0
                } else {
                    r * r * witt_mao_magnification_scaled(x * rho, r * rho)
                }
            };
            let mut result = 0.0;
            for (r, b) in bl.windows(2).zip(bf.windows(2)) {
                let slope = (b[1] - b[0]) / (r[1] - r[0]);
                let cut = x.clamp(r[0], r[1]);
                let integral = if slope == 0.0 {
                    0.0
                } else {
                    quad(&c, r[0], cut, 1e-10, 12) + quad(&c, cut, r[1], 1e-10, 12)
                };
                result += b[1] * c(r[1]) - b[0] * c(r[0]) - slope * integral;
            }
            checked(result / (2.0 * flux))
        })
        .collect()
}

fn point_source_magnification(u: f64) -> f64 {
    if u == 0.0 {
        return f64::INFINITY;
    }
    if u > 1.0 {
        let v = 1.0 / u;
        return (1.0 + 2.0 * v * v) / (1.0 + 4.0 * v * v).sqrt();
    }
    (u * u + 2.0) / (u * u.hypot(2.0))
}

fn witt_mao_magnification_scaled(l_scaled: f64, rstar_scaled: f64) -> f64 {
    let rstar_scaled2: f64 = rstar_scaled.powi(2);

    if rstar_scaled <= 0.0 {
        return point_source_magnification(l_scaled);
    }

    let l_r_diff: f64 = l_scaled - rstar_scaled;
    let l_r_sum: f64 = l_scaled + rstar_scaled;

    if l_scaled == 0.0 {
        return rstar_scaled.hypot(2.0) / rstar_scaled;
    }
    if l_r_diff == 0.0 {
        return ((2.0 / rstar_scaled)
            + ((1.0 + rstar_scaled2) / rstar_scaled2) * (2.0 * rstar_scaled.atan()))
            / f64::consts::PI;
    }

    // The elliptic terms become ill-conditioned near the limb and when the
    // disk is tiny compared with its separation. Use the exact lens-centered
    // area integral in those regions, retaining the actual separation.
    if l_r_diff.abs() < 0.01 * rstar_scaled || l_scaled > 10.0 * rstar_scaled {
        return disk_area_magnification(l_scaled / rstar_scaled, 1.0 / rstar_scaled);
    }

    let kernel1: f64 = (l_r_diff).powi(2);
    let kernel2: f64 = (4.0 + kernel1).sqrt();

    let elliptic_n: f64 = 4.0 * rstar_scaled * l_scaled / (l_r_sum).powi(2);

    let elliptic_k: f64 = (4.0 * elliptic_n).sqrt() / kernel2;
    let elliptic_m: f64 = elliptic_k.powi(2);

    let first_term: f64 = (l_r_sum * kernel2) / (2.0 * rstar_scaled2);
    let second_term: f64 =
        l_r_diff * (4.0 + (0.5 * (l_scaled.powi(2) - rstar_scaled2))) / (kernel2 * rstar_scaled2);
    let third_term: f64 =
        2.0 * kernel1 * (1.0 + rstar_scaled2) / (rstar_scaled2 * l_r_sum * kernel2);

    let ellip1: f64 = match ellip::ellipe(elliptic_m) {
        Ok(v) => v,
        Err(_e) => return f64::NAN,
    };
    let ellip2: f64 = match ellip::ellipk(elliptic_m) {
        Ok(v) => v,
        Err(_e) => return f64::NAN,
    };
    let ellip3: f64 = match ellip::ellippi(elliptic_n, elliptic_m) {
        Ok(v) => v,
        Err(_e) => return f64::NAN,
    };

    // Witt & Mao have a ±π in their equation, but we only ever want the total which
    // simplifies to the following.
    (ellip1 * first_term - ellip2 * second_term + ellip3 * third_term) / f64::consts::PI
}

// Independent form of the uniform-disk integral: in lens-centered polar
// coordinates integral s A0(s) ds = s sqrt(s^2+4 epsilon^2)/2.
fn disk_area_magnification(l: f64, e: f64) -> f64 {
    let limit = if l <= 1.0 {
        f64::consts::PI
    } else {
        (1.0 / l).asin()
    };
    let f = |theta: f64| {
        let d = (1.0 - (l * theta.sin()).powi(2)).max(0.0).sqrt();
        let center = l * theta.cos();
        let upper = center + d;
        if l <= 1.0 {
            let upper = if center < 0.0 {
                (1.0 - l * l) / (d - center)
            } else {
                upper
            };
            upper * upper.hypot(2.0 * e) / 2.0
        } else {
            let lower = center - d;
            // Rationalize the difference of the two antiderivative values.
            2.0 * center * d * (upper * upper + lower * lower + 4.0 * e * e)
                / (upper * upper.hypot(2.0 * e) + lower * lower.hypot(2.0 * e))
        }
    };
    2.0 / f64::consts::PI * quad(&f, 0.0, limit, 1e-10 * (1.0 + e), 12)
}

#[pyfunction]
pub fn witt_mao_magnification(l: Vec<f64>, re: f64, rstar: f64) -> PyResult<Vec<f64>> {
    validate(&l, re, rstar)?;
    let rstar_scaled: f64 = rstar / re;

    use rayon::prelude::*;

    l.par_iter()
        .map(|_l| {
            let l_scaled: f64 = _l * rstar_scaled;
            checked(witt_mao_magnification_scaled(l_scaled, rstar_scaled))
        })
        .collect()
}

#[pyfunction]
pub fn multi_witt_mao_magnification(
    l: Vec<f64>,
    re: Vec<f64>,
    rstar: Vec<f64>,
) -> PyResult<Vec<Vec<f64>>> {
    let mut res: Vec<Vec<f64>> = Vec::new();

    use rayon::prelude::*;

    // Prepare all combinations of re and rstar
    let combinations: Vec<(f64, f64)> = re
        .iter()
        .flat_map(|&_re| rstar.iter().map(move |&_rstar| (_re, _rstar)))
        .collect();

    // Compute in parallel
    let results: Vec<Result<Vec<f64>, PyErr>> = combinations
        .par_iter()
        .map(|&(re_val, rstar_val)| witt_mao_magnification(l.clone(), re_val, rstar_val))
        .collect();

    // Collect results, returning early on error
    for result in results {
        match result {
            Ok(v) => res.push(v),
            Err(e) => return Err(e),
        }
    }
    Ok(res)
}

/// A Python module implemented in Rust.
#[pymodule]
fn rustlens(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(witt_mao_magnification, m)?)?;
    m.add_function(wrap_pyfunction!(integrated_witt_mao_magnification, m)?)?;
    m.add_function(wrap_pyfunction!(
        integrated_flux_map_witt_mao_magnification,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(multi_witt_mao_magnification, m)?)?;

    m.add_function(wrap_pyfunction!(heyrovsky_magnification, m)?)?;
    m.add_function(wrap_pyfunction!(integrated_heyrovsky_magnification, m)?)?;
    m.add_function(wrap_pyfunction!(
        integrated_flux_map_heyrovsky_magnification,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        discrete_flux_map_heyrovsky_magnification,
        m
    )?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            actual.is_finite() && (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
            "actual={actual:.16e}, expected={expected:.16e}, relative error={:.3e}",
            (actual / expected - 1.0).abs()
        );
    }

    // Direct angular integration of the point-source formula, independent of
    // elliptic integrals, with a split resolving the closest approach.
    fn angular_reference(l: f64, r: f64, e: f64) -> f64 {
        let f = |t: f64| {
            point_source_magnification(
                ((l - r).powi(2) + 4.0 * l * r * (t / 2.0).sin().powi(2)).sqrt() / e,
            )
        };
        let split = ((l - r).abs() / (l * r).sqrt()).min(f64::consts::PI);
        (quad(&f, 0.0, split, 1e-9, 16) + quad(&f, split, f64::consts::PI, 1e-9, 16))
            / f64::consts::PI
    }

    #[test]
    fn ring_matches_point_lens_angular_integral() {
        pyo3::Python::initialize();
        for e in [0.01, 0.25, 1.0, 10.0, 100.0] {
            for (l, r) in [
                (0.2, 0.7),
                (0.7, 0.2),
                (0.5, 0.50001),
                (1e-7, 2e-7),
                (2.0, 0.4),
                (100.0, 0.8),
            ] {
                close(
                    ring_magnification(l, r, e),
                    angular_reference(l, r, e),
                    2e-8,
                );
                close(
                    ring_magnification(l, r, e),
                    ring_magnification(r, l, e),
                    1e-13,
                );
            }
        }
    }

    #[test]
    fn analytic_limits() {
        pyo3::Python::initialize();
        for e in [0.001, 0.1, 1.0, 10.0, 1000.0] {
            close(
                witt_mao_magnification_scaled(0.0, 1.0 / e),
                (1.0 + 4.0 * e * e).sqrt(),
                1e-13,
            );
            close(
                ring_magnification(0.0, 0.7, e),
                point_source_magnification(0.7 / e),
                1e-13,
            );
            close(
                ring_magnification(0.7, 0.0, e),
                point_source_magnification(0.7 / e),
                1e-13,
            );
            assert!(ring_magnification(0.5, 0.5, e).is_infinite());
            assert!(ring_magnification(0.0, 0.0, e).is_infinite());
            close(
                witt_mao_magnification_scaled(1.0 / e, 1.0 / e),
                disk_area_magnification(1.0, e),
                1e-9,
            );
        }
    }

    #[test]
    fn uniform_disks_agree_with_independent_area_integral() {
        pyo3::Python::initialize();
        for e in [1e-6, 0.001, 0.01, 0.25, 1.0, 10.0, 100.0, 1e6] {
            for l in [
                0.0, 1e-7, 0.2, 0.9, 0.999, 0.999999, 1.0, 1.000001, 1.001, 1.1, 2.0, 20.0, 1000.0,
            ] {
                let expected = disk_area_magnification(l, e);
                let disk = witt_mao_magnification_scaled(l / e, 1.0 / e);
                close(disk, expected, 2e-8);
                let ring = integrate_rings(l, e, &|_| 1.0, &[0.0, 1.0], 0.5)
                    .unwrap_or_else(|err| panic!("l={l}, e={e}: {err}"));
                close(ring, expected, 2e-8);
            }
        }
    }

    #[test]
    fn limb_darkening_methods_agree() {
        pyo3::Python::initialize();
        let positions = vec![0.0, 1e-6, 0.2, 0.8, 0.9999, 1.0, 1.0001, 1.2, 3.0, 100.0];
        for e in [0.01, 0.25, 1.0, 10.0] {
            let disk = integrated_witt_mao_magnification(positions.clone(), e, 1.0).unwrap();
            let ring = integrated_heyrovsky_magnification(positions.clone(), e, 1.0).unwrap();
            for (a, b) in disk.iter().zip(ring) {
                close(*a, b, 3e-8);
            }
            close(*disk.last().unwrap(), 1.0, 3e-4);
        }
    }

    #[test]
    fn radial_profiles_and_flux_normalization() {
        pyo3::Python::initialize();
        let positions = vec![0.0, 0.2, 0.5, 0.9999, 1.0, 1.3, 20.0];
        for (bl, bf) in [
            (vec![0.0, 1.0], vec![1.0, 1.0]),
            (vec![0.0, 0.2, 0.5, 1.0], vec![1.0, 0.8, 0.9, 0.1]),
            (vec![0.2, 0.4, 0.8], vec![0.2, 2.0, 0.7]),
        ] {
            for e in [0.1, 1.0, 10.0] {
                let a = integrated_flux_map_witt_mao_magnification(
                    positions.clone(),
                    e,
                    1.0,
                    bl.clone(),
                    bf.clone(),
                )
                .unwrap();
                let b = integrated_flux_map_heyrovsky_magnification(
                    positions.clone(),
                    e,
                    1.0,
                    bl.clone(),
                    bf.clone(),
                )
                .unwrap();
                let scaled = integrated_flux_map_witt_mao_magnification(
                    positions.clone(),
                    e,
                    1.0,
                    bl.clone(),
                    bf.iter().map(|x| 7.0 * x).collect(),
                )
                .unwrap();
                for ((a, b), c) in a.iter().zip(b).zip(scaled) {
                    close(*a, b, 3e-8);
                    close(*a, c, 1e-12);
                }
            }
        }
        let a = discrete_flux_map_heyrovsky_magnification(
            vec![3.0],
            1.0,
            1.0,
            vec![0.2, 0.7],
            vec![2.0, 3.0],
        )
        .unwrap()[0];
        close(
            a,
            (2.0 * ring_magnification(3.0, 0.2, 1.0) + 3.0 * ring_magnification(3.0, 0.7, 1.0))
                / 5.0,
            1e-13,
        );
    }

    #[test]
    fn invalid_inputs_return_errors() {
        pyo3::Python::initialize();
        assert!(witt_mao_magnification(vec![0.0], 0.0, 1.0).is_err());
        assert!(heyrovsky_magnification(-1.0, 0.2, 1.0, 1.0).is_err());
        for (bl, bf) in [
            (vec![], vec![]),
            (vec![0.0, 1.0], vec![0.0, 0.0]),
            (vec![1.0, 0.0], vec![1.0, 1.0]),
            (vec![0.0, 1.0], vec![1.0]),
            (vec![0.0, 1.0], vec![1.0, f64::NAN]),
        ] {
            assert!(integrated_flux_map_witt_mao_magnification(
                vec![0.2],
                1.0,
                1.0,
                bl.clone(),
                bf.clone()
            )
            .is_err());
            assert!(
                integrated_flux_map_heyrovsky_magnification(vec![0.2], 1.0, 1.0, bl, bf).is_err()
            );
        }
    }
}
