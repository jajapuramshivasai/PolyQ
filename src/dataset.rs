//! Dataset module for handling named polynomials.
//!
//! Provides a struct for associating a name with a polynomial object.

use crate::polynomial::Polynomial;

/// A dataset containing a name and a polynomial.
#[derive(Clone, Debug)]
pub struct Dataset {
    /// The name of the dataset.
    pub name: String,
    /// The associated polynomial.
    pub polynomial: Polynomial,
}

impl Dataset {
    /// Creates a new `Dataset` with the given name and polynomial.
    ///
    /// # Arguments
    /// * `name` - The name of the dataset.
    /// * `polynomial` - The polynomial to associate with the dataset.
    pub fn new(name: &str, polynomial: Polynomial) -> Self {
        Self {
            name: name.to_string(),
            polynomial,
        }
    }
}
