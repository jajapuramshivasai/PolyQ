use std::f64::consts::SQRT_2;
use std::fs::File;
use std::io::Write;

use num_complex::Complex64;
use numpy::PyArray1;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyModule;

type Term = (u8, Vec<usize>);

fn eval_terms(terms: &[Term], x: &[bool]) -> u8 {
    let mut val_out: u8 = 0;
    for (weight, indices) in terms {
        let mut v = true;
        for &j in indices {
            v &= x[j];
            if !v {
                break;
            }
        }
        if v {
            val_out = (val_out + (*weight % 8)) % 8;
        }
    }
    val_out
}

fn eval_terms_no_ivs(terms: &[Term], x: &[bool], n: usize) -> u8 {
    let mut val_out: u8 = 0;
    for (weight, indices) in terms {
        let mut v = true;
        for &j in indices {
            let idx = j.checked_sub(n).expect("term index underflow");
            v &= x[idx];
            if !v {
                break;
            }
        }
        if v {
            val_out = (val_out + (*weight % 8)) % 8;
        }
    }
    val_out
}

fn project_bits(index: usize, ovs: &[usize], t: usize) -> usize {
    let mut result = 0usize;
    for &bit_index in ovs {
        let shift = t - 1 - bit_index;
        let bit = (index >> shift) & 1;
        result = (result << 1) | bit;
    }
    result
}

fn accumulate_counts(ttb: &[u8], t: usize, ovs: &[usize]) -> Vec<[i64; 8]> {
    let out_size = 1usize << ovs.len();
    let mut counts = vec![[0i64; 8]; out_size];
    for (k, &t_val) in ttb.iter().enumerate() {
        let chosen_int = project_bits(k, ovs, t);
        counts[chosen_int][t_val as usize] += 1;
    }
    counts
}

fn amplitude_from_counts(counts: &[i64; 8], scale: f64) -> Complex64 {
    let tmp0 = ((counts[1] - counts[5]) as f64) / SQRT_2;
    let tmp1 = ((counts[3] - counts[7]) as f64) / SQRT_2;
    let real = (counts[0] - counts[4]) as f64 + tmp0 - tmp1;
    let imag = (counts[2] - counts[6]) as f64 + tmp0 + tmp1;
    Complex64::new(real / scale, imag / scale)
}

#[pyfunction]
fn create_poly(
    instructions: Vec<(String, Vec<usize>)>,
    n: usize,
) -> PyResult<(Vec<Term>, Vec<Vec<usize>>, usize)> {
    let mut wire_array: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut max_new_var = n;
    let mut terms: Vec<Term> = Vec::new();

    for (gate, qubits) in instructions {
        match gate.as_str() {
            "h" => {
                if qubits.is_empty() {
                    return Err(PyValueError::new_err("Hadamard gate missing target qubit"));
                }
                let target = qubits[0];
                if target >= wire_array.len() {
                    return Err(PyValueError::new_err("Qubit index out of range"));
                }
                let prev_idx = *wire_array[target]
                    .last()
                    .ok_or_else(|| PyValueError::new_err("Wire array empty"))?;
                let new_idx = max_new_var;
                wire_array[target].push(new_idx);
                max_new_var += 1;
                terms.push((4, vec![prev_idx, new_idx]));
            }
            "z" | "cz" | "ccz" => {
                let indices: PyResult<Vec<usize>> = qubits
                    .iter()
                    .map(|&q| {
                        wire_array
                            .get(q)
                            .and_then(|wire| wire.last().copied())
                            .ok_or_else(|| PyValueError::new_err("Qubit index out of range"))
                    })
                    .collect();
                terms.push((4, indices?));
            }
            "s" => {
                let indices: PyResult<Vec<usize>> = qubits
                    .iter()
                    .map(|&q| {
                        wire_array
                            .get(q)
                            .and_then(|wire| wire.last().copied())
                            .ok_or_else(|| PyValueError::new_err("Qubit index out of range"))
                    })
                    .collect();
                terms.push((2, indices?));
            }
            "t" => {
                let indices: PyResult<Vec<usize>> = qubits
                    .iter()
                    .map(|&q| {
                        wire_array
                            .get(q)
                            .and_then(|wire| wire.last().copied())
                            .ok_or_else(|| PyValueError::new_err("Qubit index out of range"))
                    })
                    .collect();
                terms.push((1, indices?));
            }
            "sdg" => {
                let indices: PyResult<Vec<usize>> = qubits
                    .iter()
                    .map(|&q| {
                        wire_array
                            .get(q)
                            .and_then(|wire| wire.last().copied())
                            .ok_or_else(|| PyValueError::new_err("Qubit index out of range"))
                    })
                    .collect();
                terms.push((6, indices?));
            }
            "tdg" => {
                let indices: PyResult<Vec<usize>> = qubits
                    .iter()
                    .map(|&q| {
                        wire_array
                            .get(q)
                            .and_then(|wire| wire.last().copied())
                            .ok_or_else(|| PyValueError::new_err("Qubit index out of range"))
                    })
                    .collect();
                terms.push((7, indices?));
            }
            other => {
                return Err(PyValueError::new_err(format!(
                    "Unsupported gate: {}",
                    other
                )));
            }
        }
    }
    Ok((terms, wire_array, max_new_var))
}

#[pyfunction]
fn get_truthtable(
    terms: Vec<Term>,
    n: usize,
    t: usize,
    initial_state: Vec<u8>,
) -> PyResult<Vec<u8>> {
    if initial_state.len() != n {
        return Err(PyValueError::new_err(
            "Initial state length must equal number of qubits",
        ));
    }
    if t < n {
        return Err(PyValueError::new_err(
            "Total variable count must be >= qubit count",
        ));
    }
    let ancilla = t - n;
    if ancilla == 0 {
        return Ok(vec![]);
    }
    let mut x = vec![false; t];
    for (i, &value) in initial_state.iter().enumerate() {
        x[i] = value != 0;
    }
    let size = 1usize << ancilla;
    let mut ttb = vec![0u8; size];
    for i in 0..size {
        for ind in 0..ancilla {
            let shift = ancilla - 1 - ind;
            x[n + ind] = ((i >> shift) & 1) != 0;
        }
        ttb[i] = eval_terms(&terms, &x);
    }
    Ok(ttb)
}

#[pyfunction]
fn get_truthtable_no_ivs(
    terms: Vec<Term>,
    n: usize,
    t: usize,
    initial_state: Vec<u8>,
) -> PyResult<Vec<u8>> {
    if initial_state.len() != n {
        return Err(PyValueError::new_err(
            "Initial state length must equal number of qubits",
        ));
    }
    if t < n {
        return Err(PyValueError::new_err(
            "Total variable count must be >= qubit count",
        ));
    }
    let ancilla = t - n;
    if ancilla == 0 {
        return Ok(vec![]);
    }
    let filtered_terms: Vec<Term> = terms
        .into_iter()
        .filter(|(_, indices)| indices.iter().all(|&idx| idx >= n))
        .collect();
    let mut x = vec![false; ancilla];
    let size = 1usize << ancilla;
    let mut ttb = vec![0u8; size];
    for i in 0..size {
        for ind in 0..ancilla {
            let shift = ancilla - 1 - ind;
            x[ind] = ((i >> shift) & 1) != 0;
        }
        ttb[i] = eval_terms_no_ivs(&filtered_terms, &x, n);
    }
    Ok(ttb)
}

#[pyfunction]
fn get_statevector(
    py: Python<'_>,
    ttb: Vec<u8>,
    n: usize,
    t: usize,
    ovs: Vec<usize>,
    starting_index: Option<usize>,
) -> PyResult<Py<PyArray1<Complex64>>> {
    let _ = starting_index.unwrap_or(0);
    if ttb.is_empty() {
        return Err(PyValueError::new_err("Truth table must not be empty"));
    }
    let scale = (SQRT_2).powi((t - n) as i32);
    let counts = accumulate_counts(&ttb, t, &ovs);
    let mut statevector = Vec::with_capacity(counts.len());
    for entry in &counts {
        statevector.push(amplitude_from_counts(entry, scale));
    }
    let array = PyArray1::from_vec_bound(py, statevector);
    Ok(array.unbind())
}

#[pyfunction]
fn get_statevector_file(
    ttb: Vec<u8>,
    n: usize,
    t: usize,
    ovs: Vec<usize>,
    path: Option<String>,
    starting_index: Option<usize>,
) -> PyResult<()> {
    let _ = starting_index.unwrap_or(0);
    if ttb.is_empty() {
        return Err(PyValueError::new_err("Truth table must not be empty"));
    }
    let scale = (SQRT_2).powi((t - n) as i32);
    let counts = accumulate_counts(&ttb, t, &ovs);
    let mut output = String::new();
    let width = ovs.len();
    for (k, entry) in counts.iter().enumerate() {
        let amp = amplitude_from_counts(entry, scale);
        let binary_k = format!("{k:0width$b}", width = width);
        output.push_str(&format!("k: {binary_k}, amp: {amp}\n"));
    }
    let filename = path.unwrap_or_else(|| "Results/demo/stvec_tmp.txt".to_string());
    let mut file = File::create(&filename).map_err(|err| {
        PyRuntimeError::new_err(format!("Failed to create file {filename}: {err}"))
    })?;
    file.write_all(output.as_bytes()).map_err(|err| {
        PyRuntimeError::new_err(format!("Failed to write file {filename}: {err}"))
    })?;
    Ok(())
}

#[pymodule]
fn polyq_backend(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(create_poly, m)?)?;
    m.add_function(wrap_pyfunction!(get_truthtable, m)?)?;
    m.add_function(wrap_pyfunction!(get_truthtable_no_ivs, m)?)?;
    m.add_function(wrap_pyfunction!(get_statevector, m)?)?;
    m.add_function(wrap_pyfunction!(get_statevector_file, m)?)?;
    Ok(())
}
