import numpy as np
import math
import scipy.sparse

fout = open( "test.metdense", "wb" )

fout.write( b"MetDense" )
fout.write( int(0).to_bytes( 4, "little" ) )  # major version
fout.write( int(0).to_bytes( 4, "little" ) )  # minor version
assert fout.tell() == 16
fout.write( int(0).to_bytes( 4, "little" ) )  # pointer to dense data
fout.write( int(0).to_bytes( 4, "little" ) )  # pointer to position information

ncells = 0
with open("/home/anders/papagei/mnt/raid/leoni/scnmt_gastrulation/scbs/data_from_leonie/column_header.txt") as fh:
    for l in fh:
        fout.write( l.encode("UTF-8") )
        ncells += 1
fout.write( b'\n' )

data_offset = fout.tell()
fout.seek( 16 )
fout.write( int(data_offset).to_bytes( 4, "little" ) )  # pointer to dense data
fout.seek( 0, 2 )  # end of file

sp = scipy.sparse.load_npz( "/home/anders/papagei/mnt/raid/leoni/scnmt_gastrulation/scbs/data_from_leonie/19.npz" )
assert ncells == sp.shape[1]
positions = []

def write_pos_to_file( fout, sp, pos ):
    res = np.empty(math.ceil(sp.shape[1] / 16), np.uint32)
    for i in range(res.shape[0]):
        a = 0
        for j in range(min(16, ncells - i * 16)):
            a <<= 2
            v = sp[pos, i * 16 + j]
            if v == 0:
                continue
            if sp[pos, i * 16 + j] == 1:
                a += 0x02
            elif sp[pos, i * 16 + j] == -1:
                a += 0x01
            else:
                raise ("Error")
        res[i] = a
    res.tofile(fout)

for pos in range(3100000):  #range(sp.shape[0]):
    if sp.indptr[pos] == sp.indptr[pos + 1]:
        continue
    positions.append(pos)
    print( pos )
    write_pos_to_file( fout, sp, pos )

print(positions)


fout.write( int(24).to_bytes( 4, "little" ) )  # pointer to dense data


fout.close()

