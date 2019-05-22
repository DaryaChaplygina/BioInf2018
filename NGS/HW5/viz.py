import pygraphviz as pgv 

if __name__ == "__main__":
    for n in ['1000', '10000', '100000']:
        for end_ in ['_no_tails', '_no_err_edges', '_compressed']:
            name = f's6_{n}{end_}'
            f = name + '.dot'
            
            g = pgv.AGraph(f, directed=False)
            g.layout(prog='dot')
            g.draw(path=f'{name}.pdf')
