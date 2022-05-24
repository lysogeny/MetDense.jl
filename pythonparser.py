import numpy as np
from numba import jit
import os
import bisect

class MetDenseFile:
    def __init__(self, filename):
        self.filename = filename
        self.size = os.path.getsize(filename)
        self.size_in_blocks = self.size//4 # in 4-blocks
        self.file = open(self.filename, "rb")
        
        self.read_intro()
        self.read_chromosomes_map()
        self.read_cell_names()
    
    ###### Basic functions
    def __len__(self): # thinks in 4-blocks
        return self.size_in_blocks
    def __getitem__(self, i): # To run binary search directly on this object, thinks in 4-blocks
        assert i < self.size_in_blocks # If index is beyond end of file
        return self.get_int(4*i)
    
    def get_int(self, position=None): # read a single int32
        if position==None:
            return int.from_bytes(self.file.read(4),"little")
        else:
            self.file.seek(position)
            return int.from_bytes(self.file.read(4),"little")
    ###### Read initial stuff
    def read_intro(self): # Read the intro part of the metdense file
        self.file.seek(8)
        self.vMajor, self.vMinor, self.data_start, self.chromosome_start, self.n_cells = [self.get_int() for i in range(5)]
        self.data_row_size = int(np.ceil(self.n_cells/16))
        self.data_row_size_bytes = 4*self.data_row_size
    
    def read_until_end(self): # Read until the file is finished
        res = b''
        cur = self.file.read(4)
        while cur:
            res += cur
            cur = self.file.read(4)
        return res
    
    def read_chromosomes_map(self): # Read the chromosome map
        self.file.seek(self.chromosome_start)
        self.n_chromosomes = self.get_int()
        self.positions_start = [self.get_int() for i in range(self.n_chromosomes)]
        self.positions_end = self.positions_start[1:]+[self.chromosome_start]
        self.data_block_size = self.positions_start[0]-self.data_start
        self.chromosome_names = self.read_until_end().split()
        self.chromosome_names = [name.decode() for name in self.chromosome_names]
        self.n_CpG = (self.positions_end[-1]-self.positions_start[0])//4
        self.n_CpG_per_chromosome = [(self.positions_end[i] - self.positions_start[i])//4 for i in range(self.n_chromosomes)]
    
    def read_cell_names(self): # Read list of cell names
        self.file.seek(28)
        self.cell_names = list()
        for i in range(self.n_cells):
            self.cell_names.append(self.file.readline()[:-1].decode())
    
    ###### Data access
    def get_position_range(self, chromosome, lower=0, upper=1e30):
        assert str(chromosome) in self.chromosome_names
        assert lower <= upper
        pos = self.chromosome_names.index(str(chromosome))
        lower_bound, upper_bound = self.positions_start[pos]//4, self.positions_end[pos]//4
        assert self[lower_bound] <= upper
        lower_index = 4*bisect.bisect_left(self, lower, lo=lower_bound, hi=upper_bound) if lower != 0 else 4*lower_bound
        if lower == upper:
            return lower_index, lower_index+4
        upper_index = 4*bisect.bisect_right(self, upper, lo=lower_bound, hi=upper_bound) if upper != 1e30 else 4*upper_bound
        return lower_index, upper_index # upper includes the index of the first CpG that is no longer included!
    
    def position_to_row_start(self, position):
        row = (position - self.positions_start[0])//4
        row_position = self.data_start + row*self.data_row_size_bytes
        return row_position
    
    #@staticmethod
    #@jit(nopython=True)
    def read_data_rows(self, file, lower, upper, row_size):
        def read_block():
            res = list()
            block = self.get_int()
            for i in range(16):
                res.append(block & 0x03)
                block >>= 2
            res.reverse()
            return res
        def read_single_row():
            row = list()
            for i in range(row_size):
                row = row + read_block()
            return row
        
        data = list()
        file.seek(lower)
        while file.tell() < upper:
            row = read_single_row()
            data.append(row)
        return data
    
    def read_data_range(self, chromosome, lower=0, upper=1e30):
        """
            - chromosome: which chromosome
            - lower: lower position on the chromosome
            - upper: upper position on the chromosome
        """
        lower_index, upper_index = self.get_position_range(chromosome, lower, upper)
        row_lower, row_upper = self.position_to_row_start(lower_index), self.position_to_row_start(upper_index)
        return np.asarray(self.read_data_rows(self.file, row_lower, row_upper, self.data_row_size))