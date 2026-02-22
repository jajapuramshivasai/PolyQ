#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
    pub coeffs: Vec<i64>,
}

impl Polynomial {
    pub fn new(mut coeffs: Vec<i64>) -> Self {
        let mut p = Polynomial { coeffs };
        p.normalize();
        p
    }

    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty() {
            0
        } else {
            self.coeffs.len() - 1
        }
    }

    pub fn eval(&self, x: i64, modulus: Option<i64>) -> i64 {
        let mut acc: i128 = 0;
        let mut pow: i128 = 1;
        let m_opt = modulus.map(|m| m as i128);
        let xx = x as i128;

        for &c in &self.coeffs {
            let term = (c as i128) * pow;
            acc += term;
            if let Some(m) = m_opt {
                acc %= m;
                pow = (pow * xx) % m;
            } else {
                pow = pow * xx;
            }
        }

        if let Some(m) = m_opt {
            let mut out = (acc % m) as i64;
            if out < 0 { out += m as i64; }
            out
        } else {
            acc as i64
        }
    }

    pub fn normalize(&mut self) {
        while let Some(&0) = self.coeffs.last() {
            self.coeffs.pop();
        }
    }
}
