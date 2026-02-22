use crate::polynomial::Polynomial;

/// Naive O(n*m) polynomial multiplication.
pub fn multiply_naive(a: &Polynomial, b: &Polynomial) -> Polynomial {
    if a.coeffs.is_empty() || b.coeffs.is_empty() {
        return Polynomial::new(vec![]);
    }
    let na = a.coeffs.len();
    let nb = b.coeffs.len();
    let mut out = vec![0i64; na + nb - 1];
    for i in 0..na {
        for j in 0..nb {
            out[i + j] = out[i + j].wrapping_add(a.coeffs[i].wrapping_mul(b.coeffs[j]));
        }
    }
    Polynomial::new(out)
}

/// Karatsuba multiplication with a small threshold.
pub fn multiply_karatsuba(a: &Polynomial, b: &Polynomial) -> Polynomial {
    let na = a.coeffs.len();
    let nb = b.coeffs.len();
    if na == 0 || nb == 0 {
        return Polynomial::new(vec![]);
    }
    // Use naive for small sizes
    if na < 32 || nb < 32 {
        return multiply_naive(a, b);
    }

    let n = std::cmp::max(na, nb);
    let half = n / 2;

    let a_low = Polynomial::new(a.coeffs.iter().cloned().take(half).collect());
    let a_high = Polynomial::new(a.coeffs.iter().cloned().skip(half).collect());
    let b_low = Polynomial::new(b.coeffs.iter().cloned().take(half).collect());
    let b_high = Polynomial::new(b.coeffs.iter().cloned().skip(half).collect());

    let z0 = multiply_karatsuba(&a_low, &b_low);
    let z2 = multiply_karatsuba(&a_high, &b_high);

    let a_sum = add(&a_low, &a_high);
    let b_sum = add(&b_low, &b_high);
    let z1 = multiply_karatsuba(&a_sum, &b_sum);

    // z1 - z2 - z0
    let mut z1_coeffs = z1.coeffs;
    ensure_len(&mut z1_coeffs, std::cmp::max(z2.coeffs.len(), z0.coeffs.len()));
    for i in 0..z2.coeffs.len() {
        z1_coeffs[i] = z1_coeffs[i].wrapping_sub(z2.coeffs[i]);
    }
    for i in 0..z0.coeffs.len() {
        z1_coeffs[i] = z1_coeffs[i].wrapping_sub(z0.coeffs[i]);
    }

    // assemble result: z0 + (z1 << half) + (z2 << (2*half))
    let mut result = vec![0i64; z0.coeffs.len().max(z1_coeffs.len() + half).max(z2.coeffs.len() + 2 * half)];
    for i in 0..z0.coeffs.len() {
        result[i] = result[i].wrapping_add(z0.coeffs[i]);
    }
    for i in 0..z1_coeffs.len() {
        result[i + half] = result[i + half].wrapping_add(z1_coeffs[i]);
    }
    for i in 0..z2.coeffs.len() {
        result[i + 2 * half] = result[i + 2 * half].wrapping_add(z2.coeffs[i]);
    }

    Polynomial::new(result)
}

fn add(a: &Polynomial, b: &Polynomial) -> Polynomial {
    let na = a.coeffs.len();
    let nb = b.coeffs.len();
    let mut out = vec![0i64; std::cmp::max(na, nb)];
    for i in 0..na { out[i] = out[i].wrapping_add(a.coeffs[i]); }
    for i in 0..nb { out[i] = out[i].wrapping_add(b.coeffs[i]); }
    Polynomial::new(out)
}

fn ensure_len(v: &mut Vec<i64>, n: usize) {
    if v.len() < n { v.resize(n, 0); }
}

/// Public helper choosing an optimized multiply implementation.
pub fn multiply(a: &Polynomial, b: &Polynomial) -> Polynomial {
    // Pick karatsuba for larger sizes
    if a.coeffs.len() >= 64 && b.coeffs.len() >= 64 {
        multiply_karatsuba(a, b)
    } else {
        multiply_naive(a, b)
    }
}
