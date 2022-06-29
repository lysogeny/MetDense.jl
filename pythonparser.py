import numpy as np
from numba import jit
import os
import bisect
import warnings
import pandas as pd
import re
    

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
        self.cells_masked = False
        
        self.__read_intro()
        self.__read_chromosomes_map()
    
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
        return self.get_int(4*i+self.nonconform_offset)
    
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
    
    def __get_lines(self, n):
        """ Read n lines from current file position.
        """
        lines = list()
        for i in range(n):
            lines.append(self.file.readline()[:-1].decode())
        return lines
    
    ###### Read initial stuff
    def __read_intro(self):
        """ Read the intro part of the metdense file.
        """
        self.file.seek(8)
        self.vMajor, self.vMinor = [self.get_int() for i in range(2)] # Version number
        self.positions_bytesize = 4 if (self.vMajor==0 and self.vMinor==0) else 8 # uint64 for positions from 0.1 onward
        self.data_start, self.chromosome_start = [self.get_int(size=self.positions_bytesize) for i in range(2)]
        if self.data_start % 4 != 0: warnings.warn("The file does not conform to the specification! Data block starts at "+str(self.data_start)+". Parser should still work.")
        self.nonconform_offset = self.data_start % 4
        self.n_cells = self.get_int()
        self.cell_names = self.__get_lines(self.n_cells)
        self.data_row_size = int(np.ceil(self.n_cells/16)) # Row size in uint32
        self.data_row_size_bytes = 4*self.data_row_size # Row size in bytes
        if self.verbose>0: print("The file follows version {}.{} and contains data for {} cells.".format(self.vMajor,self.vMinor,self.n_cells))
     
    def __read_chromosomes_map(self):
        """ Read the chromosome map.
        """
        self.file.seek(self.chromosome_start)
        self.n_chromosomes = self.get_int()
        self.positions_start = [self.get_int(size=self.positions_bytesize) for i in range(self.n_chromosomes)]
        self.positions_end = self.positions_start[1:]+[self.chromosome_start]
        self.data_block_size = self.positions_start[0]-self.data_start
        self.chromosome_names = self.__get_lines(self.n_chromosomes)
        self.n_CpG = (self.positions_end[-1]-self.positions_start[0])//4
        self.n_CpG_per_chromosome = [(self.positions_end[i] - self.positions_start[i])//4 for i in range(self.n_chromosomes)]
        if self.data_block_size/self.data_row_size_bytes != self.n_CpG: warnings.warn("The file does not conform to the specification! Data block size, number of cells and number of CpGs does not match.")
        if self.verbose>0: print("Contains Chromosomes:",self.chromosome_names)
        if self.verbose>1: print("Chromosomes have lengths:",self.n_CpG_per_chromosome,"with total length",self.n_CpG)
    
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
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
        if not lower <= upper: raise ValueError("Upper range is not larger than lower range.")
        pos = self.chromosome_names.index(str(chromosome))
        lower_bound, upper_bound = self.positions_start[pos]//4, self.positions_end[pos]//4
        lower_bound_pos, upper_bound_pos = self[lower_bound], self[upper_bound-1]
        if not lower_bound_pos <= upper: raise ValueError("Upper bound does not include anything. Chromosome "+str(chromosome)+", between "+str(lower)+" and "+str(upper)+".")
        if not upper_bound_pos >= lower: raise ValueError("Lower bound does not include anything. Chromosome "+str(chromosome)+", between "+str(lower)+" and "+str(upper)+".")
        lower_index = 4*bisect.bisect_left(self, lower, lo=lower_bound, hi=upper_bound) if lower > lower_bound_pos else 4*lower_bound
        if lower == upper:
            if not self.get_int(lower_index) == lower: raise ValueError("A single position was requested, but does not exist.")
            return lower_index, lower_index+4
        upper_index = 4*bisect.bisect_right(self, upper, lo=lower_bound, hi=upper_bound) if upper < upper_bound_pos else 4*upper_bound
        return lower_index, upper_index # upper includes the index of the first CpG that is no longer included!
    
    def __position_to_row(self, position):
        """ Turns a position in the position block into the row number in the data block.
        """
        row = (position - self.positions_start[0])//4
        return row
    
    def __position_to_row_on_chromosome(self, chromosome, position):
        """ Turns a position in the position block into the row number in the data block of the chromosome.
        """
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
        pos = self.chromosome_names.index(str(chromosome))
        row = (position - self.positions_start[pos])//4
        return row
    
    def __position_to_row_start(self, position):
        """ Turns a position in the position block into the position in the file where the corresponding row starts.
        """
        row = self.__position_to_row(position)
        row_position = self.data_start + row*self.data_row_size_bytes
        return row_position
    
    def row_to_row_start(self, row, chromosome = ""):
        """ Turns a row in the data block into the position in the file where the row starts.
        """
        start = self.data_start + row*self.data_row_size_bytes
        if chromosome:
            if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
            pos = self.chromosome_names.index(str(chromosome))
            start += (self.positions_start[pos] - self.positions_start[0])//4*self.data_row_size_bytes
        return start
    
    def __transform_position_ranges(self, chromosome, ranges):
        """ Transform a list of ranges on the chromosome, into a list of the corresponding row ranges in the data block of that chromosome.
        """
        CpG_ranges = np.asarray([self.get_position_range(chromosome, lower=min(ran[0],ran[1]), upper=max(ran[0],ran[1])) for ran in ranges])
        return self.__position_to_row_on_chromosome(chromosome, CpG_ranges), CpG_ranges.min(), CpG_ranges.max()
        
    def __transform_position_ranges_to_indices(self, chromosome, ranges):
        """ Transform a list of ranges on the chromosome, into a list of the corresponding row ranges in the data block of that chromosome.
        """
        row_ranges, pos_min, pos_max = self.__transform_position_ranges(chromosome, ranges)
        rmin, rmax = row_ranges.min(), row_ranges.max()
        return [np.linspace(l, u-1, u-l, dtype=int) for l,u in row_ranges], rmin, rmax, pos_min, pos_max
    
    def __build_load_mask(self, ranges, length):
        """ Build mask from ranges.
        """
        mask = np.full((length),False)
        for ran in ranges:
            mask[ran] = True
        return mask
    
    def __mask_ranges_transform(self, mask, ranges):
        """ Transform indices in ranges, taking into account that some parts in between will get deleted.
        """
        key = (np.linspace(0,mask.shape[0]-1,mask.shape[0])-np.invert(mask).cumsum()).astype(int)
        ranges_tf = list()
        for ran in ranges:
            ranges_tf.append(key[ran])
        return ranges_tf
    
    def get_chromosome_length(self, chromosome):
        """ Returns the CpG length of a chromosome.
        """
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
        pos = self.chromosome_names.index(str(chromosome))
        return self.n_CpG_per_chromosome[pos]
    
    def set_cell_mask(self, mask):
        """ Set cell mask to be applied to all returned data.
        """
        #if len(mask) != self.n_cells: raise ValueError("Cell mask has wrong length.")
        self.cells_masked = True
        self.cell_mask = mask
    
    #### Read Data
    def __read_data_rows(self, lower, upper, mask = None):
        """ Read a range of data rows from the data block.
            Loading is almost instantaneous, most of the time is spent on the shift and the concatenation.
            (equal parts, no idea why concatenation is so slow; allocating the memory beforehand and setting parts of the new array is even slower)
            
            lower: lower row number (inclusive)
            upper: upper row number (exclusive)
        """
        self.file.seek(0)
        data = np.reshape(np.fromfile(self.file,self.nptype,offset=lower+self.nonconform_offset ,count=(upper-lower)//4),(-1,self.data_row_size))
        if not mask is None:
            data = data[mask].copy()
        return self.__process_data_rows(data)
    
    def __process_data_rows(self, data):
        """ Read a range of data rows from the data block.
            Loading is almost instantaneous, most of the time is spent on the shift and the concatenation.
            (equal parts, no idea why concatenation is so slow; allocating the memory beforehand and setting parts of the new array is even slower)
            
            lower: lower row number (inclusive)
            upper: upper row number (exclusive)
        """
        result = list()
        for i in range(16):
            result.append(np.bitwise_and(data,0x03).astype(np.uint8))
            data = np.right_shift(data,2)
        del data # No reason to keep it in the memory
        result = [np.array_split(res,res.shape[1],axis=1) for res in result] # Reuse the variable, or just delete the old one afterwards
        result = list(map(list, zip(*result)))
        result = np.concatenate([el for ls in result for el in ls],axis=1)
        result = result[:,:self.n_cells]
        if self.cells_masked:
            result = result[:,self.cell_mask]
        return result
    
    def read_data_range(self, chromosome, lower=0, upper=1e30, mask=None):
        """ Read all data for a CpG range on a chromosome.
            
            chromosome: which chromosome
            lower: lower position on the chromosome (inclusive)
            upper: upper position on the chromosome (inclusive)
        """
        lower_index, upper_index = self.get_position_range(chromosome, lower, upper)
        row_lower, row_upper = self.__position_to_row_start(lower_index), self.__position_to_row_start(upper_index)
        return self.__read_data_rows(row_lower, row_upper, mask)
    
    def read_data_positions_from_ranges(self, chromosome, ranges):
        """ Read all data for multiple CpG ranges on a chromosome.
            
            chromosome: which chromosome
            ranges
        """
        range_indices, rmin, rmax, pos_min, pos_max = self.__transform_position_ranges_to_indices(chromosome, ranges)
        range_indices -= rmin
        mask = self.__build_load_mask(range_indices, rmax-rmin)
        
        data = self.__read_data_rows(self.row_to_row_start(rmin, chromosome), self.row_to_row_start(rmax, chromosome), mask)
        range_indices = self.__mask_ranges_transform(mask, range_indices)
        final = [data[ran] for ran in range_indices]
        
        self.file.seek(0)
        positions = np.fromfile(self.file,self.nptype,offset=pos_min+self.nonconform_offset ,count=(pos_max-pos_min)//4)[mask]
        pos = [positions[ran] for ran in range_indices]
        
        return final, pos
    
    def read_data_chromosome(self, chromosome, mask=None):
        """ Read all data on a chromosome.
        """
        return self.read_data_range(chromosome, 0, 1e30, mask)
    
    def read_data_full(self):
        """ Read the full data.
        """
        return self.__read_data_rows(self.data_start, self.positions_start[0])
    
    #### Read Positions
    def read_positions_range(self, chromosome, lower=0, upper=1e30):
        """ Read all available CpG positions within a range on a chromosome.
            
            - chromosome: which chromosome
            - lower: lower position on the chromosome (inclusive)
            - upper: upper position on the chromosome (inclusive)
        """
        lower_index, upper_index = self.get_position_range(chromosome, lower, upper)
        self.file.seek(0)
        positions = np.fromfile(self.file,self.nptype,offset=lower_index+self.nonconform_offset ,count=(upper_index-lower_index)//4)
        return positions


@jit
def tricube_mean_sample(xlist,ylist,zlist,xwidth,ywidth,xrange,yrange):
    def tricube_mean(xlist,ylist,zlist,xcent,ycent,xwidth,ywidth):
        def tricube(x,cent,width): return 70/81*(1-(2*abs(x-cent)/width)**3)**3 if 2*abs(x-cent)/width<1 else 0
        weights = list()
        for i in range(len(xlist)):
            weights.append(tricube(xlist[i],xcent,xwidth)*tricube(ylist[i],ycent,ywidth))
        result = 0
        weightsum = 0
        for i in range(len(zlist)):
            #if weights[i] != 0:
            weightsum += weights[i]
            result += weights[i]*zlist[i]
        if weightsum!=0:
            result /= weightsum
        else:
            result = np.nan
        #result /= weightsum
        return result
    xleft = np.searchsorted(xlist,xrange-xwidth/2,"left")
    xright = np.searchsorted(xlist,xrange+xwidth/2,"right")
    res = [[tricube_mean(xlist[xL:xR],ylist[xL:xR],zlist[xL:xR],x,y,xwidth,ywidth) for y in yrange] for xL, xR, x in zip(xleft,xright,xrange)]
    return np.asarray(res)


class MetDenseFileSmoother(MetDenseFile):
    """ 
    """
    def __init__(self, filename, cell_mask, pseudotime, verbose=0):
        super().__init__(filename)
        self.set_cell_mask(cell_mask)
        self.pseudotime = pseudotime
    
    def sample_smoothed(self, xlist,ylist,zlist,xwidth,ywidth,xrange,yrange, flatten = False):
        vals = tricube_mean_sample(np.asarray(xlist),np.asarray(ylist),np.asarray(zlist),xwidth,ywidth,np.asarray(xrange),np.asarray(yrange))
        x = np.repeat(np.expand_dims(xrange,1),len(yrange),axis=1)
        y = np.repeat(np.expand_dims(yrange,0),len(xrange),axis=0)
        if flatten:
            result = {"x":x.flatten(), "y":y.flatten(), "vals":vals.flatten()}
        else:
            result = {"x":x, "y":y, "vals":vals}
        return result
    
    def __get_standard_data(self, chromosome,lower,upper, repeat = True):
        dat = self.read_data_range(chromosome,lower,upper)
        pos = self.read_positions_range(chromosome,lower,upper)
        pt = self.pseudotime
        
        if repeat:
            pos = np.repeat(np.expand_dims(pos,1),dat.shape[1],axis=1)
            pt = np.repeat(np.expand_dims(pt,0),dat.shape[0],axis=0)
        
        return dat, pos, pt

    def sample_smoothed_range(self, chromosome, lower, upper, xwidth, ywidth, xrange, yrange, flatten = False):
        dat, pos, pt = self.__get_standard_data(chromosome,lower-xwidth,upper+xwidth)
        
        x_positions = pos[np.logical_or(dat==2,dat==1)]
        y_positions = pt[np.logical_or(dat==2,dat==1)]
        z_values = dat[np.logical_or(dat==2,dat==1)].astype(float)-1
        
        return self.sample_smoothed(x_positions,y_positions,z_values,xwidth,ywidth,xrange,yrange, flatten)
     
    def get_data_points(self, chromosome, lower, upper):
        dat, pos, pt = self.__get_standard_data(chromosome,lower,upper)

        methylated = {"x":pos[dat==2], "y":pt[dat==2]}
        unmethylated = {"x":pos[dat==1], "y":pt[dat==1]}
        ambiguous = {"x":pos[dat==3], "y":pt[dat==3]}
        
        return (methylated, unmethylated, ambiguous)
     
    def get_coverage(self, chromosome, lower, upper):
        dat, pos, pt = self.__get_standard_data(chromosome,lower,upper, False)
        cov = np.logical_or(dat==1,dat==2).sum(axis=1)
        
        return pos, cov
     
    def get_coverage_densesum(self, chromosome, lower_, upper_, sumover_ = 1):
        lower, upper, sumover = int(lower_), int(upper_), int(sumover_)
        pos, cov = self.get_coverage(chromosome, lower, upper)
        densecov = np.zeros(int(upper)-int(lower))
        densecov[pos-lower] = cov
        densecov = np.pad(densecov, (0,int((sumover-densecov.shape[0]%sumover)%sumover)), 'constant', constant_values=(0,0))
        densecov = densecov.reshape(-1,sumover).sum(axis=1)
        boundaries = np.asarray(list(np.arange(lower, upper,sumover))+[upper])
        
        return boundaries, densecov