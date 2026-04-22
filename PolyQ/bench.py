import time
import numpy as np
import matplotlib.pyplot as plt
from PolyQ.lib import QC

def generate_random_circuit(num_qubits, num_h, num_cz, num_rz):
    """Generate a random quantum circuit with specified parameters."""
    qc = QC(num_qubits)
    for _ in range(num_h):
        qc.h(np.random.randint(0, num_qubits))
    for _ in range(num_cz):
        q1, q2 = np.random.choice(num_qubits, 2, replace=False)
        qc.cz(q1, q2)
    for _ in range(num_rz):
        q = np.random.randint(0, num_qubits)
        angle = np.random.uniform(0, 2 * np.pi)
        qc.rz(q, angle)
    return qc

def benchmark_get_output_amplitude_0():
    """Benchmark get_output_amplitude_0 with varying parameters."""
    qubit_counts = [2, 4, 6, 8, 10, 12, 14, 16]
    h_counts = [0, 5, 10, 15, 20, 25, 30]
    cz_counts = [5, 10, 15, 20, 25]
    rz_counts = [0, 5, 10]

    results = []

    for num_qubits in qubit_counts:
        for num_h in h_counts:
            for num_cz in cz_counts:
                for num_rz in rz_counts:
                    qc = generate_random_circuit(num_qubits, num_h, num_cz, num_rz)
                    start_time = time.time()
                    qc.get_output_amplitude_0(y_val=0)
                    elapsed_time = time.time() - start_time
                    results.append((num_qubits, num_h, num_cz, num_rz, elapsed_time))

    return results, h_counts

def plot_benchmark_results(results, h_counts):
    """Plot the benchmark results and fit a polynomial."""
    qubit_counts = sorted(set(r[0] for r in results))
    avg_times = {q: [] for q in qubit_counts}

    for q in qubit_counts:
        for h in h_counts:
            filtered = [r for r in results if r[0] == q and r[1] == h]
            avg_time = np.mean([r[4] for r in filtered]) if filtered else 0
            avg_times[q].append(avg_time)

    plt.figure(figsize=(10, 6))
    for q, times in avg_times.items():
        plt.plot(h_counts, times, label=f"Qubits: {q}")

        # Fit a polynomial to the data
        poly_coeffs = np.polyfit(h_counts, times, deg=2)  # Fit a quadratic polynomial
        poly_fit = np.polyval(poly_coeffs, h_counts)

        # Print the degree of the polynomial
        print(f"Fitted polynomial for Qubits {q}: Degree = {len(poly_coeffs) - 1}")

        # Plot the polynomial fit
        plt.plot(h_counts, poly_fit, linestyle="--", label=f"Poly Fit (Qubits: {q})")

    plt.xlabel("H Gate Count")
    plt.ylabel("Execution Time (s)")
    plt.title("Scaling of get_output_amplitude_0 with Random Circuits")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    results, h_counts = benchmark_get_output_amplitude_0()
    plot_benchmark_results(results, h_counts)