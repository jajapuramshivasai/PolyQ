//TODO: split output statevector bitstring among all nodes and then combine results at the end

/*
Example: 4 qubits, 16 basis states, 4 workers
Spawn all workers with f(input,x_i),g(x_i) ,output_bitstring range specific to worker



Worker 1        || {0000 to 0011} \.  ||
Worker 2        || {0100 to 0111} /\. || final state vector
Worker 2        || {1000 to 1011} \/.  ||
Worker 2        || {1100 to 1111} /.  ||

*/
