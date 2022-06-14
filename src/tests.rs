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
        true, 1e-5, 10_000).unwrap());
}

#[test]
fn earth_to_mars_short() {
    println!("Velocities: {:?}", crate::get_velocities([1.496e11, 0.0, 0.0], [-2.28e11, 1.0e6, 0.0], 22372601.7678, 1.327e20,
        false, 1e-4, 10_000).unwrap());
}

#[test]
#[should_panic]
fn negative_tof() {
    println!("Velocities: {:?}", crate::get_velocities([1.496e11, 0.0, 0.0], [-2.28e11, 0.0, 0.0], -2.0, 1.327e20,
        false, 1e-5, 10_000).unwrap());
}

#[test]
fn ninety_degrees() {
    println!("Velocities: {:?}", crate::get_velocities([1.496e11, 0.0, 0.0], [0.0, 2.3e11, 0.0], 1e7, 1.327e20,
        false, 1e-5, 10_000).unwrap());
}

#[test]
fn bench_random() {
    use std::time::Instant;
    use rand::prelude::*;
    const NUM_TRIALS: usize = 10_000;

    let mut rng = StdRng::seed_from_u64(34585089);

    let mut r1 = vec![[0.0; 3]; NUM_TRIALS];
    let mut r2 = vec![[0.0; 3]; NUM_TRIALS];
    let mut v1s = vec![[0.0; 3]; NUM_TRIALS];
    let mut v2s = vec![[0.0; 3]; NUM_TRIALS];
    let mut tof = vec![0.0; NUM_TRIALS];
    for i in 0..NUM_TRIALS {
        r1[i][0] = (rng.gen::<f64>() - 0.5) * 1.0e11;
        r1[i][1] = (rng.gen::<f64>() - 0.5) * 1.0e11;
        r1[i][2] = (rng.gen::<f64>() - 0.5) * 1.0e11;
        r2[i][0] = (rng.gen::<f64>() - 0.5) * 1.0e11;
        r2[i][1] = (rng.gen::<f64>() - 0.5) * 1.0e11;
        r2[i][2] = (rng.gen::<f64>() - 0.5) * 1.0e11;
        tof[i] = rng.gen::<f64>() * 1.0e8;
    }
    let time_start = Instant::now();
    for i in 0..NUM_TRIALS {
        let (v1, v2) = match crate::get_velocities(r1[i], r2[i], tof[i], 1.327e20, false, 1e-7, 100){
            Ok(v) => v,
            Err(e) => panic!("{:?}, {}, {:?}, {:?}, {}", e, i, r1[i], r2[i], tof[i])
        };
        v1s[i] = v1;
        v2s[i] = v2;
    }
    let time_total = time_start.elapsed();
    //println!("{:?}", vels);
    println!("Took {:?}, {:?} us per iteration", time_total, time_total.as_micros() as f32 / NUM_TRIALS as f32);
    assert_eq!(1, 2);
}