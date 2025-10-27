#![feature(unboxed_closures)]
#![feature(fn_traits)]

use num_traits::real::Real;
use std::fmt::Display;

pub struct Polynomial<T: Real = f64> {
    coefficients: Vec<T>,
    // roots: Option<Vec<f64>>,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Degree {
    Zero,
    NonZero(usize),
}

impl Degree {
    pub fn num_degree(&self) -> Option<usize> {
        match self {
            Degree::Zero => None,
            Degree::NonZero(d) => Some(*d),
        }
    }
}

impl<T: Real> Polynomial<T> {
    pub fn new(coefficients: &[T]) -> Self {
        if coefficients.is_empty() {
            panic!("Coefficient array cannot be empty");
        }
        if coefficients.len() > 1 && coefficients[0].is_zero() {
            panic!("Leading coefficient cannot be zero for non-zero degree polynomials");
        }
        Self {
            coefficients: coefficients.to_vec(),
            // roots: None,
        }
    }

    pub fn degree(&self) -> Degree {
        if self.coefficients[0].is_zero() {
            Degree::Zero
        } else {
            Degree::NonZero(self.coefficients.len() - 1)
        }
    }

    pub fn is_zero(&self) -> bool {
        matches!(self.degree(), Degree::Zero)
    }
}

impl<T: Real + Display> Display for Polynomial<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "{}", T::zero());
        }
        let n = self.degree().num_degree().unwrap();
        let coe = &self.coefficients;
        let var = "x";
        if n == 0 {
            return write!(f, "{}", coe[0]);
        }

        let mut str = if n == 1 {
            if coe[0].abs().is_one() {
                if coe[0] > T::zero() {
                    var.to_string()
                } else {
                    format!("-{}", var)
                }
            } else {
                format!("{}{}", coe[0], var)
            }
        } else if coe[0].abs().is_one() {
            if coe[0] > T::zero() {
                format!("{}^{}", var, n)
            } else {
                format!("-{}^{}", var, n)
            }
        } else {
            format!("{}{}^{}", coe[0], var, n)
        };
        for (k, &c) in coe.iter().take(n - 1).enumerate().skip(1) {
            if c.is_zero() {
                continue;
            }
            if c > T::zero() {
                str.push('+');
            }
            if c.abs().is_one() {
                str.push_str(&format!("{}^{}", var, n - k));
            } else {
                str.push_str(&format!("{}{}^{}", c, var, n - k));
            }
        }

        if n > 1 && !coe[n - 1].is_zero() {
            if coe[n - 1] > T::zero() {
                str.push('+');
            }
            if coe[n - 1].abs().is_one() {
                str.push_str(var);
            } else {
                str.push_str(&format!("{}{}", coe[n - 1], var));
            }
        }

        if !coe[n].is_zero() {
            if coe[n] > T::zero() {
                str.push('+');
            }
            str.push_str(&format!("{}", coe[n]));
        }

        write!(f, "{}", str)
    }
}

fn evaluate_polynomial(coefficients: &[f64], x: f64) -> f64 {
    coefficients.iter().fold(0.0, |acc, &a| acc * x + a)
}

impl FnOnce<(f64,)> for Polynomial {
    type Output = f64;

    extern "rust-call" fn call_once(self, args: (f64,)) -> Self::Output {
        evaluate_polynomial(&self.coefficients, args.0)
    }
}

impl FnMut<(f64,)> for Polynomial {
    extern "rust-call" fn call_mut(&mut self, args: (f64,)) -> Self::Output {
        evaluate_polynomial(&self.coefficients, args.0)
    }
}

impl Fn<(f64,)> for Polynomial {
    extern "rust-call" fn call(&self, args: (f64,)) -> Self::Output {
        evaluate_polynomial(&self.coefficients, args.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_display() {
        let p = Polynomial::new(&[1.0, 1.0]);
        println!("{}", p);
    }
}
