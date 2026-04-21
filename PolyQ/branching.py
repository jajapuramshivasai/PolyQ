import numpy as np
import cmath
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe # Import patheffects

# --- 1. Gate Definitions ---
class Gate: pass
class H(Gate):
    def __init__(self, q): self.q = q; self.name = f"H({q})"
class Z(Gate):
    def __init__(self, q): self.q = q; self.name = f"Z({q})"
class S(Gate):
    def __init__(self, q): self.q = q; self.name = f"S({q})"
class CZ(Gate):
    def __init__(self, q1, q2): self.q1, self.q2 = q1, q2; self.name = f"CZ({q1},{q2})"
class Rz(Gate):
    def __init__(self, q, theta):
        self.q = q
        self.theta = theta
        self.name = f"Rz({q}, {theta:.2f})"

# --- 2. Pauli Decomposition Logic ---
def decompose_gate(gate):
    if isinstance(gate, (H, Z, S, CZ)):
        return [(1.0 + 0j, [gate])]
    elif isinstance(gate, Rz):
        w0 = np.cos(gate.theta / 2.0) + 0j
        w1 = -1j * np.sin(gate.theta / 2.0)
        return [(w0, []), (w1, [Z(gate.q)])]
    else:
        raise ValueError(f"Gate {gate.name} not supported.")

# --- 3. Execution Tree Node ---
class BranchNode:
    def __init__(self, accumulated_weight, gate_history=None, gate_applied=None):
        self.weight = accumulated_weight
        self.gate_history = gate_history or []
        self.gate_applied = gate_applied or "ROOT"
        self.children = []

    def add_child(self, child_node):
        self.children.append(child_node)

# --- 4. Recursive Tree Builder with Pruning ---
def build_clifford_tree(gates, current_node, threshold=1e-5):
    if not gates:
        return

    current_gate = gates[0]
    remaining_gates = gates[1:]
    branches = decompose_gate(current_gate)

    if len(branches) == 1:
        current_node.gate_history.extend(branches[0][1])
        build_clifford_tree(remaining_gates, current_node, threshold)
    else:
        for branch_weight, branch_gates in branches:
            new_weight = current_node.weight * branch_weight

            # Prune if magnitude is below threshold
            if abs(new_weight) < threshold:
                continue

            new_history = current_node.gate_history + branch_gates
            label = f"{current_gate.name} -> {'I' if not branch_gates else branch_gates[0].name}"
            child_node = BranchNode(new_weight, new_history, label)

            current_node.add_child(child_node)
            build_clifford_tree(remaining_gates, child_node, threshold)

# --- 5. NetworkX Graph Construction & Visualization ---
def build_nx_graph(root_node):
    """Converts the Node tree into a NetworkX DiGraph."""
    G = nx.DiGraph()
    id_counter = [0] # Mutable counter for unique node IDs

    def traverse(node, parent_id=None):
        current_id = id_counter[0]
        id_counter[0] += 1

        mag = abs(node.weight)
        label = f"{node.gate_applied}\n|w|={mag:.2f}"

        # Add node with magnitude as an attribute for coloring
        G.add_node(current_id, label=label, mag=mag)

        if parent_id is not None:
            G.add_edge(parent_id, current_id)

        for child in node.children:
            traverse(child, current_id)

    traverse(root_node)
    return G

def hierarchy_pos(G, root=0, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
    """Custom recursive hierarchy layout for standard trees in NetworkX."""
    pos = {root: (xcenter, vert_loc)}
    children = list(G.neighbors(root))
    if children:
        dx = width / len(children)
        nextx = xcenter - width/2 - dx/2
        for child in children:
            nextx += dx
            child_pos = hierarchy_pos(G, child, width=dx, vert_gap=vert_gap,
                                      vert_loc=vert_loc-vert_gap, xcenter=nextx)
            pos.update(child_pos)
    return pos

def visualize_tree(root_node):
    """Plots the tree, color-coding nodes by their accumulated weight."""
    G = build_nx_graph(root_node)
    pos = hierarchy_pos(G, root=0)

    # Extract attributes for plotting
    labels = nx.get_node_attributes(G, 'label')
    magnitudes = list(nx.get_node_attributes(G, 'mag').values())

    plt.figure(figsize=(12, 8))

    # Draw nodes with a colormap (e.g., 'plasma' or 'viridis')
    nodes = nx.draw_networkx_nodes(
        G, pos,
        node_size=2500,
        node_color=magnitudes,
        cmap=plt.cm.plasma,
        vmin=0.0, vmax=1.0, # Weights range from 0 to 1
        edgecolors='black'
    )

    # Draw edges
    nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=20, edge_color='gray')

    # Draw labels manually to apply patheffects
    for node, (x, y) in pos.items():
        plt.text(x, y, labels[node],
                 fontsize=9,
                 fontweight='bold',
                 color='black',
                 ha='center', va='center',
                 path_effects=[pe.withStroke(linewidth=2, foreground="white")]
                )

    # Add a colorbar to act as the legend for the weights
    cbar = plt.colorbar(nodes, pad=0.02)
    cbar.set_label('Accumulated Branch Magnitude |w|', fontsize=12)

    plt.title("Sum-over-Cliffords Execution Tree\n(Nodes colored by Path Magnitude)", fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

# ==========================================
# Execution Example
# ==========================================
if __name__ == "__main__":
    gates = [
        H(0),
        H(1),
        CZ(0, 1),
        Rz(0, np.pi/4),   # T gate
        S(1),
        Rz(1, np.pi/8),   # Smaller rotation
        H(0)
    ]

    root = BranchNode(accumulated_weight=1.0 + 0j)

    # Threshold controls how aggressively we prune low-weight branches
    build_clifford_tree(gates, root, threshold=0.2)

    print("Generating tree visualization...")
    visualize_tree(root)
