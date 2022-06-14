/// The maximum value of z for one orbit is (2 * pi)^2. Enforce this. 
const MAX_Z: f64 = 39.4784176043;

/// When z becomes too large, it should be reset to some point just inside the boundary of `MAX_Z`. If this inset is too far and z becomes too large again, it is halved each iteration.
const FIRST_RESET_DIST: f64 = 7.0;

/// Enumerates the possible ways in which the Lambert solver could fail
#[derive(Clone, Copy, Debug)]
pub enum LambertError {
    /// The algorithm could not converge in the number of iterations allowed. Raise `num_iter` or increase `eps`
    NoConvergence,
    /// NaN results were encountered. This generally occurs when assumptions on r1, r2, and tof were violated.
    NaNResult,
}

/// Find the root of a function using Newton's method. This function is specifically
/// designed for Lambert's problem and incorporates bounds on z.
pub(crate) fn find_root<Func: Fn(f64) -> f64, Slope: Fn(f64) -> f64>(func: Func, slope: Slope, eps: f64, num_iter: usize) -> Result<f64, LambertError> {
    let mut reset_last_frame = false;
    let mut reset_distance = FIRST_RESET_DIST;
    let mut z = 0.0;
    for _ in 0..num_iter {
        let f = func(z);
        if f.abs() < eps {
            // Root was found
            return Ok(z);
        }
        z -= f / slope(z);
        if z.is_nan() {
            return Err(LambertError::NaNResult);
        }
        if z > MAX_Z {
            if reset_last_frame {
                reset_distance /= 2.0;
            }
            z = MAX_Z - reset_distance;
            reset_last_frame = true;
        } else {
            reset_last_frame = false;
        }
    }
    Err(LambertError::NoConvergence)
}