#[cfg(test)]
mod tests;

const Z_DERIVATIVE_THRESHOLD: f64 = 1.0e-3;

struct LambertConvergency {
    atol: f64,
    rtol: f64,
    max_iter: usize
}

const fn factorial(num: u64) -> u64 {
    match num {
        0  => 1,
        1.. => num * factorial(num-1),
    }
}
fn dot(v1: [f64; 3], v2: [f64; 3]) -> f64 {
    v1[0]*v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
}

/// Solve Lambert's problem by computing the velocities of two points in an orbit given the points'
/// location and the time of flight between them.
/// 
/// The first element of the return value is the velocity associated with the position `r1` while
/// the second is associated with `r2`. `mu` indicates the gravitational parameter of the system, 
/// the product of Newton's Gravitational constant and the mass of the central body. It is assumed
/// that the mass of the orbiting object is negligible.
/// 
/// If the velocities are desired for trajectories that take the long path between `r1` and 
/// `r2`, then `short` should be set to `false`. This option is available for elliptical transfer
/// orbits only.
/// 
/// If `r1` or `r2` are nearly colinear (the sign of the angle between them is less than `1e-7` in 
/// absolute value), then `Err` will be returned. If the minimization fails to complete for another
/// reason, for example, if `tof` is negative, then `Err` will also be returned.
/// 
/// Even if `Ok` is returned, results could still be `NaN`.
/// 
/// `z_init` cannot be zero. This will cause `NaN`s to be returned. 
/// 
/// # Examples
/// 
/// To compute the time velocity required to leave Earth from the x axis and arrive at Mars on the -x
/// axis (a Hohmann transfer) in meters per second after 8 months (20736000 seconds): 
/// 
/// `get_velocities([1.496e11, 0.0, 0.0], [-2.28e11, 1.0e6, 0.0], 22372601.7678, 1.327e20,
///     true, -1.0, 1e-5, 10_000);`
/// 
/// A small offset was added to the positions to prevent them from being co-linear.
/// 
/// If the same function is run with `short` changed to `false`, the resulting velocities will 
/// have reversed y coordinate, except 
pub fn get_velocities(r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, short: bool, z_init: f64,
    eps: f64, max_iter: usize) -> Result<([f64; 3], [f64; 3]), ()> {
    let r1_mag = dot(r1, r1).sqrt();
    let r2_mag = dot(r2, r2).sqrt();
    let cos_Dnu = dot(r1, r2) / (r1_mag * r2_mag);
    let sin_Dnu = f64::sqrt(1.0 - cos_Dnu * cos_Dnu) * (if short {1.0} else {-1.0});
    if sin_Dnu.abs() < 1e-7 {
        return Err(());
    }

    let A = f64::sqrt(r1_mag * r2_mag) * sin_Dnu / f64::sqrt(1.0 - cos_Dnu);
    
    let best_z = match roots::find_root_newton_raphson(z_init, |z| -> f64 {
        
        let S = get_S(z);
        let C = get_C(z);
        let y = r1_mag + r2_mag - A * (1.0 - z * S) / C.sqrt(); // Equation 5.3-9
        let x = f64::sqrt(y / C); // Equation 5.3-10
        let t = (x * x * x * S + A * y.sqrt()) / mu.sqrt();// Equation 5-3.12
        // if !z.is_nan() {
        //     println!("{}, {}, {}, {}, {}, {}", z, S, C, y, x, t);
        // }
        t - tof
    }, |z| -> f64 {
        let S = get_S(z);
        let dS_dz = get_dS_dz(z);
        let C = get_C(z);
        let dC_dz = get_dC_dz(z);
        let y = r1_mag + r2_mag - A * (1.0 - z * S) / C.sqrt();
        let x = f64::sqrt(y / C);
        (x * x * x * (dS_dz - 3.0 * S * dC_dz / (2.0 * C)) 
            + A / 8.0 * (3.0 * S * y.sqrt() / C + A / x)) / mu.sqrt()
    }, &mut roots::SimpleConvergency { eps, max_iter } ) {
        Ok(z) => z,
        Err(_) => return Err(())
    };

    let S = get_S(best_z);
    let C = get_C(best_z);
    let y = r1_mag + r2_mag - A * (1.0 - best_z * S) / C.sqrt();
    let f = 1.0 - y / r1_mag; // Equation 5.3-13
    let g = A * f64::sqrt(y / mu);  // Equation 5.3-14
    let g_dot = 1.0 - y / r2_mag; // Equation 5.3-15

    let v1 = [
        (r2[0] - f * r1[0]) / g,
        (r2[1] - f * r1[1]) / g,
        (r2[2] - f * r1[2]) / g,
    ]; // Equation 5.3-16
    let v2 = [
        (g_dot * r2[0] - r1[0]) / g,
        (g_dot * r2[1] - r1[1]) / g,
        (g_dot * r2[2] - r1[2]) / g,
    ]; // Equation 5.3-17

    Ok((v1, v2))
}

/// Compute the S value defined by Bate
fn get_S(z: f64) -> f64 { // Equation 4.4-10
    if z > 0.0 {
        let z_sqrt = z.sqrt();
        (z_sqrt - z_sqrt.sin()) / (z * z_sqrt)
    } else {
        let z_sqrt = (-z).sqrt();
        (z_sqrt.sinh() - z_sqrt) / (-z * z_sqrt)
    }
}

/// Compute the C value defined by Bate
fn get_C(z: f64) -> f64 { // Equation 4.4-11
    if z > 0.0 {
        (1.0 - z.sqrt().cos()) / z
    } else {
        (1.0 - (-z).sqrt().cosh()) / z
    }
}

/// Compute the derivative of S
fn get_dS_dz(z: f64) -> f64 { // Equation 5.3-13
    if z.abs() < Z_DERIVATIVE_THRESHOLD {
        - 1.0 / factorial(5) as f64
        + 2.0 * z / factorial(7) as f64
        - 3.0 * z * z / factorial(9) as f64
        + 4.0 * z * z * z / factorial(11) as f64
        - 5.0 * z * z * z * z / factorial(13) as f64
    } else {
        if z > 0.0 {
            let sqrt_z = z.sqrt();
            - 0.5  / (z * z) * (2.0 + sqrt_z.cos()
                - 3.0 * sqrt_z.sin() / sqrt_z)
        } else {
            let sqrt_z = (-z).sqrt();
            - 0.5 / (z * z) * (2.0 + sqrt_z.cosh()
                - 3.0 * sqrt_z.sinh() / sqrt_z)
        }
    }
}

/// Compute the derivative of G
fn get_dC_dz(z: f64) -> f64 { // Equation 5.3-14
    if z.abs() < Z_DERIVATIVE_THRESHOLD {
        - 1.0 / factorial(4) as f64
        + 2.0 * z / factorial(6) as f64
        - 3.0 * z * z / factorial(8) as f64
        + 4.0 * z * z * z / factorial(10) as f64
        - 5.0 * z * z * z * z / factorial(10) as f64
    } else {
        if z > 0.0 {
            let sqrt_z = z.sqrt();
            - 0.5 / (z * z) * (2.0 - 2.0 * sqrt_z.cos()
                - sqrt_z.sin() * sqrt_z)
        } else {
            let sqrt_z = (-z).sqrt();
            - 0.5 / (z * z) * (2.0 - 2.0 * sqrt_z.cosh()
                + sqrt_z.sinh() * sqrt_z)
        }
    }
}