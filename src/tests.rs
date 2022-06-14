
#[test]
fn fact() {
    assert_eq!(crate::factorial(5), 120);
}

#[test]
fn test_dS() {
    let z = -1.1;
    let dz = 1e-5;
    let nd = (crate::get_S(z+dz) - crate::get_S(z-dz)) / (2.0 * dz);
    let ad = crate::get_dS_dz(z);
    println!("Numerical derivative: {}\tAnalytical derivative: {}", nd, ad);
    assert_eq!((nd - ad).abs() < dz, true);
}

#[test]
fn test_dC() {
    let z = -1.1;
    let dz = 1e-5;
    let nd = (crate::get_C(z+dz) - crate::get_C(z-dz)) / (2.0 * dz);
    let ad = crate::get_dC_dz(z);
    println!("Numerical derivative: {}\tAnalytical derivative: {}", nd, ad);
    assert_eq!((nd - ad).abs() < dz, true);
}

#[test]
fn earth_to_mars_long() {
    println!("Velocities: {:?}", crate::get_velocities([1.496e11, 0.0, 0.0], [-2.28e11, 1.0e6, 0.0], 22372601.7678, 1.327e20,
        true, -1.1, 1e-5, 10_000).unwrap());
}

#[test]
fn earth_to_mars_short() {
    println!("Velocities: {:?}", crate::get_velocities([1.496e11, 0.0, 0.0], [-2.28e11, 1.0e6, 0.0], 22372601.7678, 1.327e20,
        false, -1.1, 1e-5, 10_000).unwrap());
}

#[test]
#[should_panic]
fn negative_tof() {
    println!("Velocities: {:?}", crate::get_velocities([1.496e11, 0.0, 0.0], [-2.28e11, 1.0e6, 0.0], -2.0, 1.327e20,
        false, -1.1, 1e-5, 10_000).unwrap());
}