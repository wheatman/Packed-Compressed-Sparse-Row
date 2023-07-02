from pcsr import PCSR

p = PCSR(5)

p.add_edge(1,2,45)
p.add_edge(1,3,67)
p.add_edge(4,3,1)
    
p.print_graph()
print(p.sparse_matrix_vector_multiplication([2,2,2,2,2]))
print(p.pagerank([0,1,2,3,4]))
print(p.bfs(1))