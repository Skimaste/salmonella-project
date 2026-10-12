import Bio

print(Bio.__version__)

import Bio.SeqIO

class Node:
    """ Class Node to represent a vertex in the de bruijn graph """
    def __init__(self, lab):
        self.label = lab
        self.indegree = 0
        self.outdegree = 0

class Edge:
    def __init__(self, lab):
        self.label = lab

def read_real(fname, max=None):
    reads = [str(record.seq) for record in Bio.SeqIO.parse(fname,"fasta")]
    if max:
        return reads[:max]
    return reads

def read_reads(fname):
    """ Read short reads in FASTA format. It is assumed that one line in the input file correspond to one read. """
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    reads = []

    for line in lines:
        if line[0] != '>':
            reads = reads + [line.rstrip()]

    return reads

def construct_graph(reads, k):
    """ Construct de bruijn graph from sets of short reads with k length word"""
    edges = dict()
    vertices = dict()

    for read in reads:
        i = 0
        while i+k < len(read):
            v1 = read[i:i+k]
            v2 = read[i+1:i+k+1]
            if v1 in edges.keys():
                vertices[v1].outdegree += 1
                edges[v1] += [Edge(v2)]
            else:
                vertices[v1] = Node(v1)
                vertices[v1].outdegree += 1
                edges[v1] = [Edge(v2)]
            if v2 in edges.keys():
                vertices[v2].indegree += 1
            else:
                vertices[v2] = Node(v2)
                vertices[v2].indegree += 1
                edges[v2] = []
            i += 1

    return (vertices, edges)

def output_contigs(g):
    """Return Eulerian contigs while preserving every edge in each component.

    Components with multiple start/end imbalances are split into several
    contigs by adding temporary edges, running Hierholzer's algorithm, and
    cutting the resulting circuit at those temporary edges.
    """
    V, E = g
    adjacency = {
        label: [edge.label for edge in E.get(label, [])]
        for label in V
    }

    # An Eulerian path is checked per weakly connected component. The reverse
    # links make the component search include edges in either direction.
    undirected = {label: set() for label in V}
    for source, destinations in adjacency.items():
        for destination in destinations:
            undirected[source].add(destination)
            undirected[destination].add(source)

    remaining = {
        label for label, destinations in adjacency.items() if destinations
    }
    contigs = []

    while remaining:
        component = set()
        pending = [remaining.pop()]
        while pending:
            current = pending.pop()
            if current in component:
                continue
            component.add(current)
            for neighbour in undirected[current]:
                if neighbour not in component:
                    remaining.discard(neighbour)
                    pending.append(neighbour)

        starts = []
        ends = []
        for label in component:
            imbalance = V[label].outdegree - V[label].indegree
            if imbalance > 0:
                starts.extend([label] * imbalance)
            elif imbalance < 0:
                ends.extend([label] * -imbalance)

        if len(starts) != len(ends):
            raise ValueError("Component has unbalanced total degree")

        # Add virtual end -> start edges so the component has an Eulerian
        # circuit. Each virtual edge will become a contig boundary later.
        component_edges = {
            label: [(destination, False) for destination in adjacency[label]]
            for label in component
        }
        for end, start in zip(ends, starts):
            component_edges[end].append((start, True))

        start = next(label for label in component if component_edges[label])

        # Hierholzer traversal. Entries record the edge used to enter a node;
        # reversing the dead-end order produces the complete Eulerian circuit.
        stack = [(start, None)]
        circuit = []
        while stack:
            current = stack[-1][0]
            if component_edges[current]:
                destination, is_virtual = component_edges[current].pop()
                stack.append((destination, is_virtual))
            else:
                circuit.append(stack.pop())
        circuit.reverse()

        expected_edges = sum(len(adjacency[label]) for label in component)
        if len(circuit) != expected_edges + len(starts) + 1:
            raise ValueError("Component could not be traversed completely")

        cycle_vertices = [label for label, _ in circuit[:-1]]
        cycle_edges = [is_virtual for _, is_virtual in circuit[1:]]
        virtual_positions = [
            index for index, is_virtual in enumerate(cycle_edges)
            if is_virtual
        ]
        if not virtual_positions:
            path = cycle_vertices
            path.append(path[0])
            contigs.append(path[0] + "".join(label[-1] for label in path[1:]))
            continue

        for position in virtual_positions:
            path = [cycle_vertices[(position + 1) % len(cycle_vertices)]]
            index = (position + 1) % len(cycle_edges)
            while not cycle_edges[index]:
                path.append(cycle_vertices[(index + 1) % len(cycle_vertices)])
                index = (index + 1) % len(cycle_edges)
            if len(path) > 1:
                contigs.append(path[0] + "".join(label[-1] for label in path[1:]))

    return contigs

def print_graph(g):
    """ Print the information in the graph to be (somewhat) presentable """
    V = g[0]
    E = g[1]
    for k in V.keys():
        print(f"name: ", V[k].label, ". indegree: ", V[k].indegree, ". outdegree: ", V[k].outdegree)
        print("Edges: ")
        for e in E[k]:
            print(e.label)
        print

if __name__ == "__main__":
    fname = "./data/reference-data/salmonella-enterica.reads.fna"

    reads = read_real(fname)

    print(len(reads[0]))
    print(len(reads))
    

    k = 45

    print(k)

    g = construct_graph(reads, k)

    contigs = output_contigs(g)
    contigs.sort(key=len, reverse=True)

    for i, contig in enumerate(contigs):
        if i >= 10:
            break
        print(len(contig))
        print(contig)