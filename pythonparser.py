import numpy as np
from numba import jit
import os
import bisect

class MetDenseFile:
    """ Class to load data from a .metdense file.
        The most relevant method is read_data_range, which can load the methylation data from within a range of positions on a chromosome.
        
        filename: path to the file
        verbose: verbosity of the initialization, from 0-2
    """
    def __init__(self, filename, verbose=0):
        self.filename = filename
        self.size = os.path.getsize(filename)
        self.size_in_blocks = self.size//4 # in 4-blocks
        self.file = open(self.filename, "rb")
        self.verbose = verbose
        if self.verbose>1: print("Load file",self.filename)
        self.endianess = "little" # Should also work with "big", but untested!
        self.nptype = "<u4" if self.endianess == "little" else ">u4"
        
        self.read_intro()
        self.read_chromosomes_map()
    
    ###### Basic functions
    def __len__(self):
        """ Length, thinks in 4-blocks.
        """
        return self.size_in_blocks
    
    def __getitem__(self, i):
        """ To run binary search directly on this object, thinks in 4-blocks.
            
            i: address of the value
        """
        if not i < self.size_in_blocks: raise IndexError("Index is beyond end of file.")
        return self.get_int(4*i)
    
    def get_int(self, position=None, size=4):
        """ Read a single number from file.
            
            position: defaults to current, otherwise byte position in file
            size: number of bytes to read into int, by default 4 (uint32)
        """ # 
        if position==None:
            return int.from_bytes(self.file.read(size),self.endianess)
        else:
            self.file.seek(position)
            return int.from_bytes(self.file.read(size),self.endianess)
    
    def get_lines(self, n):
        """ Read n lines from current file position.
        """
        lines = list()
        for i in range(n):
            lines.append(self.file.readline()[:-1].decode())
        return lines
    
    ###### Read initial stuff
    def read_intro(self):
        """ Read the intro part of the metdense file.
        """
        self.file.seek(8)
        self.vMajor, self.vMinor = [self.get_int() for i in range(2)] # Version number
        self.positions_bytesize = 4 if (self.vMajor==0 and self.vMinor==0) else 8 # uint64 for positions from 0.1 onward
        self.data_start, self.chromosome_start = [self.get_int(size=self.positions_bytesize) for i in range(2)]
        self.n_cells = self.get_int()
        self.cell_names = self.get_lines(self.n_cells)
        self.data_row_size = int(np.ceil(self.n_cells/16)) # Row size in uint32
        self.data_row_size_bytes = 4*self.data_row_size # Row size in bytes
        if self.verbose>0: print("The file follows version {}.{} and contains data for {} cells.".format(self.vMajor,self.vMinor,self.n_cells))
     
    def read_chromosomes_map(self):
        """ Read the chromosome map.
        """
        self.file.seek(self.chromosome_start)
        self.n_chromosomes = self.get_int()
        self.positions_start = [self.get_int(size=self.positions_bytesize) for i in range(self.n_chromosomes)]
        self.positions_end = self.positions_start[1:]+[self.chromosome_start]
        self.data_block_size = self.positions_start[0]-self.data_start
        self.chromosome_names = self.get_lines(self.n_chromosomes)
        self.n_CpG = (self.positions_end[-1]-self.positions_start[0])//4
        self.n_CpG_per_chromosome = [(self.positions_end[i] - self.positions_start[i])//4 for i in range(self.n_chromosomes)]
        if self.verbose>0: print("Contains Chromosomes:",self.chromosome_names)
        if self.verbose>1: print("Chromosomes have lengths:",self.n_CpG_per_chromosome)
    
    ###### Data access
    #### Helpers
    def get_position_range(self, chromosome, lower=0, upper=1e30):
        """ Returns the position range in the positions block for a chromosome CpG range.
            
            chromosome: name of the chromosome
            lower: lower range on the chromosome (inclusive)
            upper: upper range on the chromosome (inclusive)
            
            lower_index: index of the first position in the range
            upper_index: index of the first position above the range (or on the next chromosome)
        """
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome does not exist.")
        if not lower <= upper: raise ValueError("Upper range is not larger than lower range.")
        pos = self.chromosome_names.index(str(chromosome))
        lower_bound, upper_bound = self.positions_start[pos]//4, self.positions_end[pos]//4
        lower_bound_pos, upper_bound_pos = self[lower_bound], self[upper_bound]
        if not lower_bound_pos <= upper: raise ValueError("Upper bound does not include anything.")
        if not upper_bound_pos >= lower: raise ValueError("Lower bound does not include anything.")
        lower_index = 4*bisect.bisect_left(self, lower, lo=lower_bound, hi=upper_bound) if lower > lower_bound_pos else 4*lower_bound
        if lower == upper:
            if not self.get_int(lower_index) == lower: raise ValueError("A single position was requested, but does not exist.")
            return lower_index, lower_index+4
        upper_index = 4*bisect.bisect_right(self, upper, lo=lower_bound, hi=upper_bound) if upper < upper_bound_pos else 4*upper_bound
        return lower_index, upper_index # upper includes the index of the first CpG that is no longer included!
    
    def position_to_row_start(self, position):
        """ Turns a position in the position block into the position in the file where the corresponding row starts.
        """
        row = (position - self.positions_start[0])//4
        row_position = self.data_start + row*self.data_row_size_bytes
        return row_position
    
    #### Read Data
    def read_data_rows(self, lower, upper):
        """ Read a range of data rows from the data block.
            
            lower: lower row number (inclusive)
            upper: upper row number (exclusive)
        """
        self.file.seek(0)
        data = np.reshape(np.fromfile(self.file,self.nptype,offset=lower ,count=(upper-lower)//4),(-1,self.data_row_size))
        result = list()
        for i in range(16):
            result.append(np.bitwise_and(data,0x03).astype(np.uint8))
            data = np.right_shift(data,2)
        del data # No reason to keep it in the memory
        result = [np.array_split(res,res.shape[1],axis=1) for res in result] # Reuse the variable, or just delete the old one afterwards
        result = list(map(list, zip(*result)))
        result = np.concatenate([el for ls in result for el in ls],axis=1)
        result = result[:,:self.n_cells]
        return result
    
    def read_data_range(self, chromosome, lower=0, upper=1e30):
        """ Read all data for a CpG range on a chromosome.
            
            chromosome: which chromosome
            lower: lower position on the chromosome (inclusive)
            upper: upper position on the chromosome (inclusive)
        """
        lower_index, upper_index = self.get_position_range(chromosome, lower, upper)
        row_lower, row_upper = self.position_to_row_start(lower_index), self.position_to_row_start(upper_index)
        return self.read_data_rows(row_lower, row_upper)
    
    def read_data_full(self):
        """ Read the full data.
        """
        return self.read_data_rows(self.data_start, self.self.positions_start[0])
    
    #### Read Positions
    def read_positions_range(self, chromosome, lower=0, upper=1e30):
        """ Read all available CpG positions within a range on a chromosome.
            
            - chromosome: which chromosome
            - lower: lower position on the chromosome (inclusive)
            - upper: upper position on the chromosome (inclusive)
        """
        lower_index, upper_index = self.get_position_range(chromosome, lower, upper)
        self.file.seek(0)
        positions = np.fromfile(self.file,self.nptype,offset=lower_index ,count=(upper_index-lower_index)//4)
        return positions