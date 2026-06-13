//! GF(256) arithmetic using the polynomial x^8 + x^4 + x^3 + x + 1
//! (the same polynomial used in AES, standard for Reed-Solomon).

use serde::{Serialize, Deserialize};

/// A GF(256) element backed by a u8.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GF256(pub u8);

// Generator polynomial for log/exp tables: x^8 + x^4 + x^3 + x + 1 = 0x11B
const PRIMITIVE: u16 = 0x11D;

// Build log and exp tables at compile time.
const fn build_tables() -> ([u8; 256], [u8; 255]) {
    let mut exp = [0u8; 255];
    let mut log = [0u8; 256];
    let mut val: u16 = 1;
    let mut i = 0;
    while i < 255 {
        exp[i] = val as u8;
        log[val as usize] = i as u8;
        val <<= 1;
        if val >= 256 {
            val ^= PRIMITIVE;
        }
        i += 1;
    }
    (log, exp)
}

static TABLES: ([u8; 256], [u8; 255]) = build_tables();
static LOG_TABLE: &[u8; 256] = &TABLES.0;
static EXP_TABLE: &[u8; 255] = &TABLES.1;

impl GF256 {
    pub const ZERO: GF256 = GF256(0);
    pub const ONE: GF256 = GF256(1);

    #[inline]
    pub fn new(v: u8) -> Self {
        GF256(v)
    }

    #[inline]
    pub fn add(self, other: GF256) -> GF256 {
        GF256(self.0 ^ other.0)
    }

    #[inline]
    pub fn sub(self, other: GF256) -> GF256 {
        // Same as add in GF(2^m)
        self.add(other)
    }

    #[inline]
    pub fn mul(self, other: GF256) -> GF256 {
        if self.0 == 0 || other.0 == 0 {
            return GF256::ZERO;
        }
        let log_sum = (LOG_TABLE[self.0 as usize] as usize + LOG_TABLE[other.0 as usize] as usize) % 255;
        GF256(EXP_TABLE[log_sum])
    }

    #[inline]
    pub fn div(self, other: GF256) -> GF256 {
        if other.0 == 0 {
            panic!("division by zero in GF(256)");
        }
        if self.0 == 0 {
            return GF256::ZERO;
        }
        let log_diff = (LOG_TABLE[self.0 as usize] as usize + 255 - LOG_TABLE[other.0 as usize] as usize) % 255;
        GF256(EXP_TABLE[log_diff])
    }

    #[inline]
    pub fn inv(self) -> GF256 {
        if self.0 == 0 {
            panic!("inverse of zero in GF(256)");
        }
        let log_val = LOG_TABLE[self.0 as usize] as usize;
        GF256(EXP_TABLE[(255 - log_val) % 255])
    }

    #[inline]
    pub fn pow(self, exp: u32) -> GF256 {
        if exp == 0 {
            return GF256::ONE;
        }
        if self.0 == 0 {
            return GF256::ZERO;
        }
        let log_val = LOG_TABLE[self.0 as usize] as u64;
        let result = (log_val * exp as u64) % 255;
        GF256(EXP_TABLE[result as usize])
    }

    /// Evaluate a polynomial (highest degree first) at this element.
    pub fn eval_poly(coeffs: &[GF256], x: GF256) -> GF256 {
        let mut result = GF256::ZERO;
        for &c in coeffs {
            result = result.mul(x).add(c);
        }
        result
    }

    /// Polynomial multiplication in GF(256).
    pub fn poly_mul(a: &[GF256], b: &[GF256]) -> Vec<GF256> {
        let mut result = vec![GF256::ZERO; a.len() + b.len() - 1];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                result[i + j] = result[i + j].add(ai.mul(bj));
            }
        }
        result
    }

    /// Polynomial division in GF(256). Returns (quotient, remainder).
    pub fn poly_div(dividend: &[GF256], divisor: &[GF256]) -> (Vec<GF256>, Vec<GF256>) {
        let mut rem = dividend.to_vec();
        let d_len = divisor.len();
        if d_len == 0 {
            panic!("division by zero polynomial");
        }
        let q_len = if rem.len() >= d_len { rem.len() - d_len + 1 } else { 0 };
        let mut quotient = vec![GF256::ZERO; q_len];
        let inv_lead = divisor[0].inv();

        for i in 0..q_len {
            if rem[i] == GF256::ZERO {
                continue;
            }
            let coeff = rem[i].mul(inv_lead);
            quotient[i] = coeff;
            for j in 0..d_len {
                rem[i + j] = rem[i + j].sub(coeff.mul(divisor[j]));
            }
        }

        // Trim trailing zeros from remainder
        let rem_trimmed: Vec<GF256> = rem.iter()
            .rev()
            .skip_while(|&&x| x == GF256::ZERO)
            .cloned()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect();

        (quotient, if rem_trimmed.is_empty() { vec![GF256::ZERO] } else { rem_trimmed })
    }

    /// Sum of a slice of GF(256) elements.
    pub fn sum(vals: &[GF256]) -> GF256 {
        vals.iter().fold(GF256::ZERO, |a, &b| a.add(b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_sub_same() {
        let a = GF256(0x57);
        let b = GF256(0x83);
        assert_eq!(a.add(b), a.sub(b));
    }

    #[test]
    fn test_mul_identity() {
        let a = GF256(0x53);
        assert_eq!(a.mul(GF256::ONE), a);
    }

    #[test]
    fn test_mul_zero() {
        let a = GF256(0x42);
        assert_eq!(a.mul(GF256::ZERO), GF256::ZERO);
    }

    #[test]
    fn test_div_inverse() {
        let a = GF256(0x73);
        assert_eq!(a.mul(a.inv()), GF256::ONE);
    }

    #[test]
    fn test_pow() {
        let a = GF256(0x02);
        // 2^8 = 1 in GF(256) (generator element property not guaranteed, but multiplicative group order is 255)
        assert_eq!(a.pow(255), GF256::ONE);
    }
}
