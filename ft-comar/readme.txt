This software is free for non-commercial use.

FT-COMAR: fault tolerant three-dimensional structure reconstruction
from protein contact maps
M. Vassura, L. Margara, P. Di Lena, F. Medri, P. Fariselli and
R. Casadio

bioinformatics.cs.unibo.it/FT-COMAR



Quick start

The FT-COMAR executable needs as parameters the
name of the file of the input contact map, the threshold and a binary
value.
The file containing the input contact map must have the number of
residues in the first line, then the values of the contact map (0=no
contact, 1=contact, -1=unknown). The format is the same of the web
server Matrix except that the size is in the first line of the contact map
file. You can find examples at the web server site. The threshold is the value
used to compute the contact map given, in the examples it is 12. The
binary value has to be 1 if you want to apply the simple common
neighbors filter to the input contact map, 0 otherwise. So to
reconstruct the structure of the file 1poh.cmap_sized without using filter you
run the program with the following command line:

FT-COMAR 1poh.cmap_sized 12 0

which gives as output

chain 1 er 0 time 0.02 erS 0

together with the 4 files: nat.cmap natf.cmap rec.cmap rec.pdb

- chain 1 means that the distance between residues
  which are consecutive in the sequence is less than 4 and that the
  distance between each residue and all the others is greater than 3.5 
- er 0 means that the number of differences between the input and the
  reconstructed contact map is 0
- time is the computation time
- erS is the number of differences between the input contact map and the
  reconstructed one considering also filtered areas (in this case filter
  is not used so er=erS)
- nat.cmap is the native contact map with 3 instead of -1
- natf.cmap is the native contact map filtered with 3 instead of -1 and
  3 in filtered positions
- rec.cmap is the contact map of the reconstructed structure
- rec.pdb is the coordinate file of the reconstructed structure
  formatted as the ATOM section of a pdb file


EVAcon format

To input contact maps in EVAcon format you must convert it to the
Matrix format explained above. One way is to use the EVAconvert
utility here included. It needs as parameters the length of the
protein and the EVAcon file. As example if the file pred.eva is the
predicted contact map of a protein of size 263 in EVAcon format then
the command line

EVAconvert 263 pred.eva >pred.cmap

writes a corresponding contact map in the Matrix format. Such contact
map is computed considering at most 5*263-100 contact probabilities to
be contacts (in decreasing order). All other entries of the matrix
with sequence separation >4 are considered to be 0.

As example to compute the 3D reconstruction of prediction pred.eva
extimating contact probabilities at threshold 8, and relative to a protein of length 263:

echo 263 >pred.cmap
EVAconvert 263 pred.eva >>pred.cmap
FT-COMAR pred.cmap 8 0 myrec.pdb

so that the reconstructed structure is saved in the file myrec.pdb


Contact information

bioinformatics.cs.unibo.it/FT-COMAR

For any question please contact:
vassura@cs.unibo.it