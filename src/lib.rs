use core::f64;
use ellip;
use interp::interp;
use interp::InterpMode;
use numdiff::central_difference::sderivative;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use quadrature::integrate;
use std::iter::zip;

const RADIUS_LIM: f64 = 1e-5;

const LD_COEFF: f64 = 0.6;
// The integral of ld_linear; no need to compute this every time
const LD_LINEAR_INT: f64 = 2.5132741228717936;

fn ld_linear(r: f64) -> f64 {
    return 1.0 - LD_COEFF * (1.0 - (1.0 - r.powi(2)).sqrt());
}

#[pyfunction]
pub fn integrated_heyrovsky_magnification(l: Vec<f64>, re: f64, rstar: f64) -> PyResult<Vec<f64>> {
    return _integrated_heyrovsky_magnification(l, re, rstar, &ld_linear, LD_LINEAR_INT);
}

#[pyfunction]
pub fn discrete_flux_map_heyrovsky_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    bl: Vec<f64>,
    bf: Vec<f64>,
) -> PyResult<Vec<f64>> {
    assert_eq!(bl.len(), bf.len());

    let mut res = Vec::new();
    for _l in l.iter() {
        let mut sum: f64 = 0.0;
        for (_bl, _bf) in zip(bl.iter(), bf.iter()) {
            let a = heyrovsky_magnification(_l.clone(), _bl.clone(), re, rstar)?;
            sum += a * _bf;
        }
        res.push(sum);
    }
    return Ok(res);
}

#[pyfunction]
pub fn integrated_flux_map_heyrovsky_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    bl: Vec<f64>,
    bf: Vec<f64>,
) -> PyResult<Vec<f64>> {
    assert_eq!(bl.len(), bf.len());

    let flux_map = |r: f64| -> f64 {
        return interp(&bl, &bf, r, &InterpMode::Constant(0.0));
    };

    let b_int = integrate(
        |r: f64| -> f64 { 2.0 * f64::consts::PI * r * flux_map(r) },
        0.0,
        1.0,
        1e-16,
    )
    .integral;

    return _integrated_heyrovsky_magnification(l, re, rstar, &flux_map, b_int);
}

fn _integrated_heyrovsky_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    b: &dyn Fn(f64) -> f64,
    b_int: f64,
) -> PyResult<Vec<f64>> {
    /*
    Known issues: This implementation follows the brightness integral from Witt & Mao
                  instead of Heyrovsky (see known issues with the Heyrovsky
                  implementation below). But, this produces irratic variations in
                  magnification rather than a smooth magnification profile.

     */
    let mut res = Vec::new();
    for _l in l.iter() {
        let mag_deriv = |r: f64| -> f64 {
            if r < 0.0 {
                return 0.0;
            }
            return sderivative(
                &|x: f64| -> f64 {
                    return match heyrovsky_magnification(_l.clone(), x, re, rstar) {
                        Ok(v) => v,
                        Err(_e) => f64::NAN,
                    };
                },
                r,
                None,
            );
        };
        let mag_int = integrate(
            |r: f64| -> f64 {
                let a = match heyrovsky_magnification(_l.clone(), r, re, rstar) {
                    Ok(v) => v,
                    Err(_e) => f64::NAN,
                };
                return 2.0 * f64::consts::PI * r * (a + r / 2.0 * mag_deriv(r)) * b(r);
            },
            0.0,
            1.0,
            1e-16,
        )
        .integral;

        res.push(mag_int / b_int);
    }
    return Ok(res);
}

// fn _integrated_heyrovsky_magnification(
//     l: Vec<f64>,
//     re: f64,
//     rstar: f64,
//     b: &dyn Fn(f64) -> f64,
// ) -> PyResult<Vec<f64>> {
//     /*
//     Known issues: Does not converge on 1 for large l values.
//     */
//     let mut res = Vec::new();
//     for _l in l.iter() {
//         let igrand = |r: f64| -> f64 {
//             let a = match heyrovsky_magnification(_l.clone(), r, re, rstar) {
//                 Ok(v) => v,
//                 Err(_e) => return f64::NAN,
//             };
//             return b(r) * a * r;
//         };
//         let int = integrate(&igrand, 0.0, 1.0, 1e-16).integral;
//         res.push(int);
//     }
//     return Ok(res);
// }

#[pyfunction]
pub fn heyrovsky_magnification(l: f64, r: f64, re: f64, rstar: f64) -> PyResult<f64> {
    /*
    Known issues: l=0, and r=0 is not handled correctly. Implementing this as per the paper results
                  in a division by zero. Also magnifications very close to the centre are very
                  large.
     */
    let epsilon: f64 = re / rstar;
    let epsilon2: f64 = epsilon.powi(2);
    let l_r_diff: f64 = r - l;

    if l_r_diff.abs() < RADIUS_LIM {
        // if l_r_diff < 0.0 {
        //     return Ok(0.0);
        // }
        let l2: f64 = l.powi(2);

        let term1: f64 = 2.0 * epsilon / l * (1.0 - l_r_diff / (2.0 * l));
        let term2: f64 = ((8.0 * epsilon * l) / (l_r_diff.abs() * (l2 + epsilon2).sqrt())).ln();
        let term3: f64 = 4.0 * (l / epsilon).atan();
        let term4: f64 = (epsilon * (2.0 * l2 + epsilon2)) / (l2 * (l2 + epsilon2)) * l_r_diff;
        return Ok((term1 * term2 + term3 + term4) / (2.0 * f64::consts::PI));
    }
    let l_r_sum: f64 = l + r;
    let l_r_diff2: f64 = l_r_diff.powi(2);
    let epsilon2_4: f64 = 4.0 * epsilon2;

    let term1: f64 = 4.0 / (l_r_sum * (l_r_diff2 + epsilon2_4).sqrt());
    let elliptic_k: f64 = ((4.0 * epsilon) / l_r_sum) * ((l * r) / (l_r_diff2 + epsilon2_4)).sqrt();
    let elliptic_m: f64 = elliptic_k.powi(2);
    let elliptic_n: f64 = 4.0 * l * r / (l_r_sum).powi(2);
    let ellip2: f64 = match ellip::ellipk(elliptic_m) {
        Ok(v) => v,
        Err(e) => return Err(PyRuntimeError::new_err(e)),
    };
    let ellip3: f64 = match ellip::ellippi(elliptic_n, elliptic_m) {
        Ok(v) => v,
        Err(e) => return Err(PyRuntimeError::new_err(e)),
    };
    return Ok((term1 * (2.0 * epsilon2 * ellip2 + l_r_diff2 * ellip3)) / (2.0 * f64::consts::PI));
}

#[pyfunction]
pub fn integrated_witt_mao_magnification(l: Vec<f64>, re: f64, rstar: f64) -> PyResult<Vec<f64>> {
    return _integrated_witt_mao_magnification(l, re, rstar, &ld_linear, LD_LINEAR_INT);
}

#[pyfunction]
pub fn integrated_flux_map_witt_mao_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    bl: Vec<f64>,
    bf: Vec<f64>,
) -> PyResult<Vec<f64>> {
    assert_eq!(bl.len(), bf.len());

    let flux_map = |r: f64| -> f64 {
        return interp(&bl, &bf, r, &InterpMode::Constant(0.0));
    };

    let b_int = integrate(
        |r: f64| -> f64 { 2.0 * f64::consts::PI * r * flux_map(r) },
        0.0,
        1.0,
        1e-16,
    )
    .integral;

    return _integrated_witt_mao_magnification(l, re, rstar, &flux_map, b_int);
}

fn _integrated_witt_mao_magnification(
    l: Vec<f64>,
    re: f64,
    rstar: f64,
    b: &dyn Fn(f64) -> f64,
    b_int: f64,
) -> PyResult<Vec<f64>> {
    let mut res = Vec::new();
    for umag in witt_mao_magnification(l, re, rstar)? {
        let radial_witt_mao_magnification = |r: f64| -> f64 {
            if r < 0.0 {
                return 0.0;
            }
            return umag;
        };
        let mag_deriv = |r: f64| -> f64 {
            if r < 0.0 {
                return 0.0;
            }
            return sderivative(
                &|x: f64| -> f64 { radial_witt_mao_magnification(x) },
                r,
                None,
            );
        };
        let mag_int = integrate(
            |r: f64| -> f64 {
                2.0 * f64::consts::PI
                    * r
                    * (radial_witt_mao_magnification(r) + r / 2.0 * mag_deriv(r))
                    * b(r)
            },
            0.0,
            1.0,
            1e-16,
        )
        .integral;
        res.push(mag_int / b_int);
    }
    return Ok(res);
}

#[pyfunction]
pub fn witt_mao_magnification(l: Vec<f64>, re: f64, rstar: f64) -> PyResult<Vec<f64>> {
    let rstar_scaled: f64 = rstar / re;
    let rstar_scaled2: f64 = rstar_scaled.powi(2);

    let mut res: Vec<f64> = Vec::new();

    for _l in l.iter() {
        let l_scaled: f64 = _l * rstar_scaled;

        let l_r_diff: f64 = l_scaled - rstar_scaled;
        let l_r_sum: f64 = l_scaled + rstar_scaled;

        if l_r_diff.abs() < RADIUS_LIM {
            res.push(
                ((2.0 / rstar_scaled)
                    + ((1.0 + rstar_scaled2) / rstar_scaled2)
                        * ((f64::consts::PI / 2.0)
                            + ((rstar_scaled2 - 1.0) / (rstar_scaled2 + 1.0)).asin()))
                    / f64::consts::PI,
            );
            continue;
        }

        let kernel1: f64 = (l_r_diff).powi(2);
        let kernel2: f64 = (4.0 + kernel1).sqrt();

        let elliptic_n: f64 = 4.0 * rstar_scaled * l_scaled / (l_r_sum).powi(2);

        let elliptic_k: f64 = (4.0 * elliptic_n).sqrt() / kernel2;
        let elliptic_m: f64 = elliptic_k.powi(2);

        let first_term: f64 = (l_r_sum * kernel2) / (2.0 * rstar_scaled2);
        let second_term: f64 = l_r_diff * (4.0 + (0.5 * (l_scaled.powi(2) - rstar_scaled2)))
            / (kernel2 * rstar_scaled2);
        let third_term: f64 =
            2.0 * kernel1 * (1.0 + rstar_scaled2) / (rstar_scaled2 * l_r_sum * kernel2);

        let ellip1: f64 = match ellip::ellipe(elliptic_m) {
            Ok(v) => v,
            Err(e) => return Err(PyRuntimeError::new_err(e)),
        };
        let ellip2: f64 = match ellip::ellipk(elliptic_m) {
            Ok(v) => v,
            Err(e) => return Err(PyRuntimeError::new_err(e)),
        };
        let ellip3: f64 = match ellip::ellippi(elliptic_n, elliptic_m) {
            Ok(v) => v,
            Err(e) => return Err(PyRuntimeError::new_err(e)),
        };

        let kernel3: f64 = ellip1 * first_term - ellip2 * second_term + ellip3 * third_term;

        let pos: f64 = (kernel3 + f64::consts::PI) / (2.0 * f64::consts::PI);
        let neg: f64 = (kernel3 - f64::consts::PI) / (2.0 * f64::consts::PI);
        res.push(pos + neg);
    }
    return Ok(res);
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
