# Qiskit vs PolyQ Benchmark Summary

Benchmark timings collected from the available `bench/` artifacts for the `test_q29_h10_t*_s10_z40_cz10.qasm` circuits.

| Circuit | polyq Amplitude sim time | polyq sv sim time | PolyQ parallel sv sim time | qiskit aer sv sim time |
|---|---|---|---|---|
| `test_q29_h10_t1_s10_z40_cz10.qasm` | 4.792 µs | 23.534594417 s | 13.695729375 s | 32.756588936 s |
| `test_q29_h10_t2_s10_z40_cz10.qasm` | 30 µs | 23.1752425 s | 8.191665875 s | 30.866231680 s |
| `test_q29_h10_t3_s10_z40_cz10.qasm` | 29.208 µs | 25.064738083 s | 9.096594292 s | 35.778780937 s |
| `test_q29_h10_t4_s10_z40_cz10.qasm` | 22.5 µs | 23.012990709 s | 8.080117584 s | 36.785545826 s |
| `test_q29_h10_t5_s10_z40_cz10.qasm` | n/a | n/a | n/a | 36.132961988 s |

> Notes: PolyQ benchmark values were extracted from `bench/polyq-core-result.txt`. The JSON `siminfo` files provide the recorded aggregate simulation times for each circuit on macOS arm64 with 8 CPU cores.