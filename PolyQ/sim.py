import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from qiskit import QuantumCircuit
from engine import QC as DicksonEngine # The Dickson simulator from engine.py

class BranchNode:
    """Represents a node in the universal execution tree."""
    def __init__(self, weight, gate_history=None, gate_applied="ROOT"):
        self.weight = weight
        self.gate_history = gate_history or []
        self.gate_applied = gate_applied
        self.children = []

    def add_child(self, child_node):
        self.children.append(child_node)

class UniversalQC:
    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.num_qubits = circuit.num_qubits
        self.root = BranchNode(weight=1.0 + 0j)

    def _decompose_instr(self, instr):
        """Decomposes Qiskit gates into (weight, Clifford_gate_list)."""
        name = instr.operation.name.lower()
        indices = [self.circuit.find_bit(q).index for q in instr.qubits]
        
        # Debug: Print gate being decomposed
        print(f"Decomposing gate: {name}, indices: {indices}")
        
        # Clifford logic (Single branch)
        if name in ['h', 'z', 's', 'sdg', 'cz', 'id']:
            return [(1.0 + 0j, [(name, indices)])]
        
        # Non-Clifford: Rz(theta) = cos(theta/2)I - i*sin(theta/2)Z
        elif name == 'rz':
            theta = instr.operation.params[0]
            return [
                (np.cos(theta/2), [('id', indices)]),
                (-1j * np.sin(theta/2), [('z', indices)])
            ]
        else:
            raise ValueError(f"Gate {name} is not supported.")

    def build_tree(self, threshold=1e-5):
        """Recursively builds the execution tree with threshold pruning."""
        def _build(gate_idx, current_node):
            if gate_idx == len(self.circuit.data):
                return

            instr = self.circuit.data[gate_idx]
            branches = self._decompose_instr(instr)

            # If only one branch (Clifford), collapse history into current node
            if len(branches) == 1:
                current_node.gate_history.extend(branches[0][1])
                _build(gate_idx + 1, current_node)
            else:
                for b_weight, b_gates in branches:
                    new_weight = current_node.weight * b_weight
                    if abs(new_weight) < threshold: continue # Pruning

                    label = f"{instr.operation.name} -> {b_gates[0][0]}"
                    child = BranchNode(new_weight, current_node.gate_history + b_gates, label)
                    current_node.add_child(child)
                    _build(gate_idx + 1, child)

        self.root = BranchNode(weight=1.0 + 0j)
        _build(0, self.root)

    def get_top_n_paths(self, n):
        """Extracts the N paths with the highest absolute weight."""
        all_paths = []
        def _collect(node):
            if not node.children:
                all_paths.append((node.weight, node.gate_history))
                return
            for child in node.children: _collect(child)
        
        _collect(self.root)
        return sorted(all_paths, key=lambda x: abs(x[0]), reverse=True)[:n]

    def get_amplitude(self, y_val, x_val=0, top_n=None):
        """Sums amplitudes across paths using the Dickson Engine."""
        paths = self.get_top_n_paths(top_n) if top_n else self.get_top_n_paths(len(self.root.children))
        total_amp = 0j
        
        for weight, history in paths:
            qc_cliff = QuantumCircuit(self.num_qubits)
            for g_name, indices in history:
                if g_name != 'id': getattr(qc_cliff, g_name)(*indices)
            
            engine = DicksonEngine(qc_cliff)
            engine.compile()
            total_amp += weight * engine.get_amplitude(y_val, x_val)
        return total_amp

    def visualize(self):
        """Visualizes the tree using the NetworkX logic from branching.py."""
        G = nx.DiGraph()
        id_gen = iter(range(100000))

        def traverse(node, p_id=None):
            c_id = next(id_gen)
            mag = abs(node.weight)
            G.add_node(c_id, label=f"{node.gate_applied}\n|w|={mag:.2f}", mag=mag)
            if p_id is not None: G.add_edge(p_id, c_id)
            for child in node.children: traverse(child, c_id)

        traverse(self.root)
        pos = self._hierarchy_pos(G, 0)
        
        plt.figure(figsize=(10, 6))
        mags = [G.nodes[n]['mag'] for n in G.nodes]
        nodes = nx.draw_networkx_nodes(G, pos, node_size=2000, node_color=mags, cmap=plt.cm.plasma)
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrowsize=15)
        
        labels = nx.get_node_attributes(G, 'label')
        for node, (x, y) in pos.items():
            plt.text(x, y, labels[node], fontsize=8, ha='center', va='center',
                     path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        
        plt.colorbar(nodes).set_label('|Weight|')
        plt.title("Universal Simulator Execution Tree")
        plt.axis('off')
        plt.show()

    def _hierarchy_pos(self, G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
        """Layout helper for tree visualization."""
        pos = {root: (xcenter, vert_loc)}
        children = list(G.neighbors(root))
        if children:
            dx = width / len(children)
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos.update(self._hierarchy_pos(G, child, width=dx, vert_gap=vert_gap, 
                                               vert_loc=vert_loc-vert_gap, xcenter=nextx))
        return pos