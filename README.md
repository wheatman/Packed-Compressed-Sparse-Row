# Packed-Compressed-Sparse-Row

Dynamic data structure for sparse graphs.

compile with "g++ -std=c++11" PCSR.cpp

```
  // initialize the structure with how many nodes you want it to start with
  PCSR pcsr = PCSR(10);

  // add some edges
  for (int i = 0; i < 5; i++) {
    pcsr.add_edge(i, i, 1);
  }
  // update the values of some edges

  for (int i = 0; i < 5; i++) {
    pcsr.add_edge_update(i, i, 2);
  }

  // print out the graph
  pcsr.print_graph();
```

For more information see https://ieeexplore.ieee.org/abstract/document/8547566

Please cite as:
```
@inproceedings{wheatman2018packed,
  title={Packed Compressed Sparse Row: A Dynamic Graph Representation},
  author={Wheatman, Brian and Xu, Helen},
  booktitle={2018 IEEE High Performance extreme Computing Conference (HPEC)},
  pages={1--7},
  year={2018},
  organization={IEEE}
}
```

# Subsequent Papers and Projects
This work was continued and made parallel in the paper [A Parallel Packed Memory Array to Store Dynamic Graphs](https://epubs.siam.org/doi/abs/10.1137/1.9781611976472.3)

The Parallel PMA was later used in the larger Terrace system.  More details can be found in [Terrace: A Hierarchical Graph Container for Skewed Dynamic Graphs](https://dl.acm.org/doi/abs/10.1145/3448016.3457313). The code for terrace which includes the Parallel PMA can be found at https://github.com/PASSIONLab/terrace.

Many ideas from this work also went into the creation of [Streaming Sparse Graphs using Efficient Dynamic Sets](https://ieeexplore.ieee.org/abstract/document/9671836) for which the code can be found at https://github.com/wheatman/SSTGraph.

All of these systems are parallel and much faster than the original PCSR system, but are more complex.  

PCSR is not being updated, but I will review any pull requests. 
