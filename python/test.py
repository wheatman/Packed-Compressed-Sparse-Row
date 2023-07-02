from pcsr import PCSR

p = PCSR(10)

for i in range(0,10):
    p.add_edge(i,i,i)
    
p.print_graph()