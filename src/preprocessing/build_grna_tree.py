import gzip
import dill
import random
from tqdm import tqdm

class Tree():
    __slots__ = ['A', 'T', 'C', 'G']
    def __init__(self):
        self.A = None
        self.T = None
        self.C = None
        self.G = None

    def __repr__(self):
        return f"Tree(A={self.A}, T={self.T}, C={self.C}, G={self.G})"
    
    def add_spacer(self, spacer):
        if len(spacer) != 20:
            raise ValueError("spacer length must be 20")
        spacer = spacer.upper()
        curr_node = self
        for nuc in spacer:
            if nuc not in "ATCG":
                raise ValueError("spacer must only contain A, C, T, G")
            if getattr(curr_node, nuc) is None:
                setattr(curr_node, nuc, Tree())
            curr_node = getattr(curr_node, nuc)
    
    def check_spacer(self, spacer):
        if len(spacer) != 20:
            raise ValueError("spacer length must be 20")
        spacer = spacer.upper()
        curr_node = self
        for nuc in spacer:
            if nuc not in "ATCG":
                raise ValueError("spacer must only contain A, C, T, G")
            curr_node = getattr(curr_node, nuc)
            if curr_node is None:
                return False
        return True



def set_to_tree():
    with gzip.open('./src/genome_files/dup_guides.pkl.gz',  'rb') as handle:
        bad_spacers = dill.load(handle)
    print(f"non-unique, bad spacers: {len(bad_spacers)}")

    head = Tree()
    for spacer in tqdm(bad_spacers):
        head.add_spacer(spacer)
    
    # create 100 random spacer sequences
    for _ in range(100000):
        spacer = ''.join([random.choice('ACGT') for _ in range(20)])
        if spacer in bad_spacers:
            print(spacer)
        if head.check_spacer(spacer) != (spacer in bad_spacers):
            print(f'{spacer} not in tree')

    print("exporting non-unique spacers TREE with pickle")
    with gzip.open('./src/genome_files/dup_guide_tree.pkl',  'wb') as handle:
        dill.dump(head, handle)    

set_to_tree()

# IT SEEMS LIKE THIS IS NO BETTER AT LIMITING MEMORY/STORAGE THAN SETS