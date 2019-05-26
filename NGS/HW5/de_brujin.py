import numpy as np
import warnings
import copy
from Bio import SeqIO
from sortedcontainers import SortedList, SortedDict

warnings.simplefilter('ignore')

class Vertex:
    def __init__(self, seq):
        self.e_in = []
        self.e_out = []
        self.seq = seq
    
    def __eq__(self, other):
        return self.seq == other.seq
    
    def __lt__(self, other):
        return self.seq < other.seq 
    def __le__(self, other):
        return self.seq <= other.seq
    def __ne__(self, other):
        return self.seq != other.seq
    def __gt__(self, other):
        return self.seq > other.seq
    def __ge__(self, other):
        return self.seq >= other.seq
    
    
class Edge:
    def __init__(self, start, end, seq, cov=1):
        self.start = start
        self.end = end
        self.seq = seq
        self.cov = cov
        
class DeBrujinGraph:
    def __init__(self, k : int):
        self.k = k
        self.v_list = []
        self.v_indexes = SortedDict()
        self.e_list = []
        self.v_removed = []
        self.e_removed = []
        
    def build(self, file : str, fformat = "fasta"):
        with open(file, "r") as f:
            for read in SeqIO.parse(f, fformat):
                self.__add_read(str(read.seq))
                self.__add_read(str(read.seq.reverse_complement()))
                
    def build_from_str(self, seq):
        self.__add_read(seq)
        
    def compress(self):
        for i in range(len(self.v_list)):
            if i not in self.v_removed:
                v = self.v_list[i]
                if len(v.e_in) == len(v.e_out) == 1: 
                    self.v_list[self.e_list[v.e_out[0]].end].e_in.remove(v.e_out[0])
                    
                    self.__add_edge(v.e_in[0], 
                                    self.e_list[v.e_in[0]].start, 
                                    self.e_list[v.e_out[0]].end,
                                    self.e_list[v.e_in[0]].seq + self.e_list[v.e_out[0]].seq[self.k:])
                    
                    self.e_removed.append(v.e_out[0])
                    self.v_removed.append(i)
                
    def update_edge_coverage(self, kmers_cov):
        for i in range(len(self.e_list)):
            if i not in self.e_removed:
                e = self.e_list[i] 
                
                self.e_list[i].cov = 0
                for kmer in self.__kmers(self.k + 1, e.seq):
                    self.e_list[i].cov += kmers_cov[kmer]
                self.e_list[i].cov /= (len(e.seq) - self.k)
                
    def remove_bad_edges(self):
        len_thr = self.k * 2 + 5
        
        for i in range(len(self.e_list)):
            if i not in self.e_removed:
                e = self.e_list[i]
                if len(e.seq) < len_thr and e.cov < 30:
                    self.e_removed.append(i)
                    
        self.__update_v_set()
                    
    def remove_tails(self):
        len_thr = self.k * 2
        
        for v in self.__get_vertexes():
            tail_idx = -1
            other_cov = 0
            if len(v.e_out) == 0 and len(v.e_in) == 1:
                tail_idx = v.e_in[0]
                if len(self.v_list[self.e_list[tail_idx].start].e_in) + \
                  len(self.v_list[self.e_list[tail_idx].start].e_out) < 3:
                    tail_idx = -1
                else:
                    for i in self.v_list[self.e_list[tail_idx].start].e_in:
                        if i != tail_idx:
                            other_cov += self.e_list[i].cov
                    for i in self.v_list[self.e_list[tail_idx].start].e_out:
                        if i != tail_idx:
                            other_cov += self.e_list[i].cov
                    other_cov /= len(self.v_list[self.e_list[tail_idx].start].e_in) + \
                        len(self.v_list[self.e_list[tail_idx].start].e_out)
            elif len(v.e_in) == 0 and len(v.e_out) == 1:
                tail_idx = v.e_out[0]
                if len(self.v_list[self.e_list[tail_idx].end].e_in) + \
                  len(self.v_list[self.e_list[tail_idx].end].e_out) < 3:
                    tail_idx = -1
                else:
                    for i in self.v_list[self.e_list[tail_idx].end].e_in:
                        if i != tail_idx:
                            other_cov += self.e_list[i].cov
                    for i in self.v_list[self.e_list[tail_idx].end].e_out:
                        if i != tail_idx:
                            other_cov += self.e_list[i].cov
                    other_cov /= len(self.v_list[self.e_list[tail_idx].end].e_in) + \
                        len(self.v_list[self.e_list[tail_idx].end].e_out)
                    
            if tail_idx != -1 and len(self.e_list[tail_idx].seq) < len_thr \
              and self.e_list[tail_idx].cov < other_cov:
                e = self.e_list[tail_idx]
                self.v_list[e.end].e_in.remove(tail_idx)
                self.v_list[e.start].e_out.remove(tail_idx)
                self.e_removed.append(tail_idx)
                if len(self.v_list[e.end].e_in) == len(self.v_list[e.end].e_out) == 0:
                    self.v_removed.append(e.end)
                if len(self.v_list[e.start].e_in) == len(self.v_list[e.start].e_out) == 0:
                    self.v_removed.append(e.start)
                
    def to_dot(self, out : str):
        f = open(out, "w")
        f.write(f"graph deBrujin_kmers{self.k}"+"{ \n")
        for v in self.__get_vertexes():
            f.write(f"\t{v.seq};\n")
        for e in self.__get_edges():
            f.write(f"\t{self.v_list[e.start].seq} -- {self.v_list[e.end].seq} [label=coverage_{int(e.cov)}];\n")
        f.write("}")
        
    def to_fasta(self, out: str):
        n = 1
        f = open(out, "w")
        for e in self.__get_edges(): 
            f.write(f"> edge{n}; coverage: {e.cov}\n")
            f.write(f"{e.seq}\n")
            n += 1
                
    def __add_read(self, seq):
        for kmer in self.__kmers(self.k + 1, seq):
            v1 = kmer[: self.k]
            v2 = kmer[1 :]
            
            idx1, v1_exist = self.__add_vertex(v1)
            idx2, v2_exist = self.__add_vertex(v2)  

            e_exist = False
            if v1_exist and v2_exist:
                e1 = self.v_list[idx1].e_out
                e2 = self.v_list[idx2].e_in
                intersections = set(e1) & set(e2)
                if len(intersections) == 1:
                    for inters in intersections:
                        self.e_list[inters].cov += 1
                    e_exist = True

            if not e_exist:
                self.__add_edge(len(self.e_list), idx1, idx2, kmer)
                
    def __update_v_set(self):
        for i in range(len(self.v_list)):
            if i not in self.v_removed:
                to_remove = []
                for e_idx in self.v_list[i].e_in:
                    if e_idx in self.e_removed:
                        to_remove.append(e_idx)
                for e in to_remove:
                    self.v_list[i].e_in.remove(e)
                
                to_remove = []
                for e_idx in self.v_list[i].e_out:
                    if e_idx in self.e_removed:
                        to_remove.append(e_idx)
                for e in to_remove:
                    self.v_list[i].e_out.remove(e)
                
                if len(self.v_list[i].e_in) == len(self.v_list[i].e_out) == 0:
                    self.v_removed.append(i)
        
    @staticmethod
    def __kmers(k : int, seq : str):
        i = 0
        while i <= len(seq) - k:
            yield seq[i : i + k]
            i += 1
            
    def __add_vertex(self, v : str):
        idx = self.v_indexes.get(v)
        flag = True
        if idx is None:
            idx = len(self.v_list)
            flag = False
            self.v_list.append(Vertex(v))
            self.v_indexes[v] = idx
            
        return idx, flag
    
    def __add_edge(self, e_idx, start, end, seq, cov=1):
        if e_idx == len(self.e_list):
            self.e_list.append(Edge(start, end, seq, cov))
            self.v_list[start].e_out.append(e_idx)
        else:
            self.e_list[e_idx] = Edge(start, end, seq, cov)
            
        self.v_list[end].e_in.append(e_idx)
        
    def __str__(self):
        s = ""
        for e in self.__get_edges():
            s += f'{self.v_list[e.start].seq}\t{self.v_list[e.end].seq}\tcov:{e.cov}\t{e.seq}\n'
        return s
    
    def __get_edges(self):
        for i in range(len(self.e_list)):
            if i not in self.e_removed:
                yield self.e_list[i]
                
    def __get_vertexes(self):
        for i in range(len(self.v_list)):
            if i not in self.v_removed:
                yield self.v_list[i]
                
def get_kmers_cov(file, fformat, k):
    res = SortedDict()
    with open(file, "r") as f:
        for read in SeqIO.parse(file, fformat):
            for i in range(len(read.seq) - k + 1):
                kmer = str(read.seq[i : i + k])
                if kmer in res.keys():
                    res[kmer] += 1
                else:
                    res[kmer] = 1
            for i in range(len(read.seq) - k + 1):
                kmer = str(read.seq.reverse_complement()[i : i + k])
                if kmer in res.keys():
                    res[kmer] += 1
                else:
                    res[kmer] = 1
    return res
                
if __name__ == "__main__":
    f = input()
    fformat = input()
    name_ = input()
    k = int(input())
    
    kmers_cov = get_kmers_cov(f, fformat, k + 1)
    
    g = DeBrujinGraph(k)
    #g.build_from_str("ACGTCCGTAA")
    g.build(f, fformat)
    #g.to_dot(name_ + "_uncompressed.dot")
    g.compress()
    g.update_edge_coverage(kmers_cov)
    #g.to_dot(name_ + "_compressed.dot")
    #g.to_fasta(name_ + "_compressed.fasta")
    
    """
    g1 = copy.deepcopy(g)
    prev_edges = len(g1.e_list)
    curr_edges = len(g1.e_list) - 1
    while prev_edges != curr_edges:
        g1.remove_tails()
        g1.compress()
        g1.update_edge_coverage(kmers_cov)
        prev_edges = curr_edges
        curr_edges = len(g1.e_list)
        
    g1.to_dot(name_ + "_no_tails.dot")
    g1.to_fasta(name_ + "_no_tails.fasta")
    g1 = DeBrujinGraph(k)    
    """
    
    g.remove_bad_edges()
    g.compress()
    g.update_edge_coverage(kmers_cov)
    g.to_dot(name_ + "_no_err_edges.dot")
    g.to_fasta(name_ + "_no_err_edges.fasta")

