# Packed Compressed Sparse Row for Python

The present README offers comprehensive instructions for utilizing the Python-based implementation of PCSR.

## Compiling PCSR.cpp into the .so file

Execute the following command to obtain the shared object file, which can be imported into any Python program.

```bash
c++ -O3 -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) pcsr.cpp -o pcsr$(python3-config --extension-suffix)
```

```python3
from pcsr import PCSR
```

## Using the PCSR Python Module

### Initialising a PCSR Graph Object

The following code-snippet will initialise a PCSR graph object with 4 nodes

```python
from pcsr import PCSR

graph = PCSR(4)
```

### Adding edges to the graph

To add edges into the PCSR graph, specify the source node, destination node, and the corresponding edge value using the ```PCSR.add_edge(src, dst, val)``` method.

```python
from pcsr import PCSR

graph = PCSR(4)

graph.add_edge(2,3,4)
graph.add_edge(0,1,10)
graph.add_edge(1,3,9)
```

### Printing the graph

To visualize the PCSR graph in the form of an adjacency matrix, you can utilize the ```PCSR.print_graph()``` method.

```python
from pcsr import PCSR

graph = PCSR(4)

graph.add_edge(2,3,4)
graph.add_edge(0,1,10)
graph.add_edge(1,3,9)

graph.print_graph()
```
**Output**
```
   0   1   2   3
0      10
1              9
2              4
3
```

### Adding a new node

To introduce a new node to the PCSR graph, employ the ```PCSR.add_node()``` method. This operation adds a single node identified by ```node_id```, where the ```node_id``` value is one greater than the maximum ```node_id``` present in the current graph.

```python
from pcsr import PCSR

graph = PCSR(4)

graph.add_edge(2,3,4)
graph.add_edge(0,1,10)
graph.add_edge(1,3,9)

graph.add_node()

graph.print_graph()
```
**Output**
```
   0   1   2   3   4
0      10
1              9
2              4
3
4
```
Node with ```node_id``` = 4 was added to the graph.

### Updating an edge value

The ```PCSR.add_edge_update(src, dst, val)``` method allows for the updating of an edge value. If the specified edge does not currently exist in the graph, this method will add a new edge with the provided source, destination, and value to the graph.

```python
from pcsr import PCSR

graph = PCSR(4)

graph.add_edge(2,3,4)
graph.add_edge(0,1,10)
graph.add_edge(1,3,9)

graph.add_node()

graph.add_edge_update(0,1,5)    # this edge exists
graph.add_edge_update(0,3,15)   # this edge doesn't exists

graph.print_graph()
```
**Output**
```
   0   1   2   3   4
0      5       15
1              9
2              4
3
4
```

## Get the edge list of the graph

The edge list representation of the graph can be obtained by utilizing the ```PCSR.get_edges()``` method. This method retrieves the edges of the graph in the ```[src, dst, val]``` structured format.

```python
from pcsr import PCSR

graph = PCSR(4)

graph.add_edge(2,3,4)
graph.add_edge(0,1,10)
graph.add_edge(1,3,9)

graph.add_node()

graph.add_edge_update(0,1,5)    
graph.add_edge_update(0,3,15)   

edge_list = graph.get_edges()
print(edge_list)
```
**Output**
```
[(0, 1, 5), (0, 3, 15), (1, 3, 9), (2, 3, 4)]
```

### Number of nodes and size of the graph

To obtain the number of nodes in the graph, employ the ```PCSR.get_n()``` method. Similarly, the size of the graph can be acquired using the ```PCSR.get_size()``` method.

```python

from pcsr import PCSR

graph = PCSR(4)

graph.add_edge(2,3,4)
graph.add_edge(0,1,10)
graph.add_edge(1,3,9)

graph.add_node()

graph.add_edge_update(0,1,5)    
graph.add_edge_update(0,3,15)   

print(f'Number of nodes: {graph.get_n()}')
print(f'Size of the graph: {graph.get_size()}')
```
**Output**
```
Number of nodes: 5
Size of the graph: 256
```