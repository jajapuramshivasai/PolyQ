//WIP
//todo: seperate poly into Z8+z4+z2 which can be sent to their respective engine
use std::collections::HashMap;
use crate::Z8::engine::PhasePolynomialZ8;

pub fn optimize_polynomial(poly: &mut PhasePolynomialZ8) {
    // 1. Combine Duplicate Terms
    // Use a HashMap to store the canonical representation of variable indices
    let mut term_map: HashMap<Vec<usize>, i32> = HashMap::new();
    
    for (weight, mut idxs) in poly.terms.drain(..) {
        idxs.sort_unstable(); // Ensure [1, 2] and [2, 1] are treated as the same term
        let entry = term_map.entry(idxs).or_insert(0);
        *entry = (*entry + weight) % 8;
    }

    // 2. Filter out terms with zero weight (mod 8)
    poly.terms = term_map
        .into_iter()
        .filter(|(_, weight)| *weight != 0)
        .map(|(idxs, weight)| (weight, idxs))
        .collect();

    // 3. Variable Pruning 
    // Identify internal variables (index >= num_qubits) that no longer appear in any terms
    let mut active_vars = std::collections::HashSet::new();
    for (_, idxs) in &poly.terms {
        for &idx in idxs {
            active_vars.insert(idx);
        }
    }
    
    // Update wire_array to reflect that some intermediate variables might have vanished
    // Note: In a production engine, you would re-index variables here to keep the range contiguous.
}

/*
Before running the heavy simulation loop, we can perform Symbolic Reduction to eliminate redundant or internal variables. This step prunes the polynomial by identifying variables that can be simplified analytically, drastically reducing the effective 2 
free_bits
  search space.

1. Linear Variable Elimination

If a variable xi appears only once in the entire polynomial in a linear term (e.g., 4xi), and it is not an output variable, 
its contribution to the sum over all paths is always zero unless that specific phase is balanced. More importantly, if x i
​	
appears in a term like 4xi⋅xj, 
it acts as a constraint that can often nullify other terms.

2. Term Consolidation

Many terms in the polynomial might be duplicates or can be combined using Z 
8
​	
  arithmetic (e.g., 2x1+2x1=4x1(mod8)).
After pruning variables using the methods above, the indices in poly.terms will have "gaps" (e.g., variables 0, 1, 4, 7).
Optimization: Re-map all active variables to a contiguous range [0,new_t].

Polynomial optimizations

Implement Number Theoretic Transform (NTT) multiplication with NTT-friendly moduli.
Provide Karatsuba (already added) and FFT-based multiplication fallbacks.

*/