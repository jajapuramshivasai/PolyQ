/// Small ring arithmetic helpers used by the PolyQ port.
pub fn add_mod(a: i64, b: i64, modulus: i64) -> i64 {
    let m = modulus;
    (((a % m) + (b % m)) % m + m) % m
}

pub fn sub_mod(a: i64, b: i64, modulus: i64) -> i64 {
    add_mod(a, -b, modulus)
}

pub fn mul_mod(a: i64, b: i64, modulus: i64) -> i64 {
    let m = modulus as i128;
    let res = ((a as i128) * (b as i128)) % m;
    let mut out = res as i64;
    if out < 0 { out += modulus; }
    out
}
