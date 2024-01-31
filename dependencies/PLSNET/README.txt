***********************************************************************

SOURCES


Here is a brief description of the MATLAB functions contained in the
 directory plsnet.m : 
Returns a matrix where the element (i,j) is the weight of the edge directed from the ith gene to the jth gene.

In the MATLAB prompt 

plsnet_single.m:
Returns a vector where the ith element is the weight of the edge directed from the ith gene to the target gene.
In the MATLAB prompt 

get_link_list.m: 
Writes the ranked list of inferred regulatory links.

**********************************************************************DEMO
D is a 100*100 matrix
V = plsnet(D,1:100,5,30,1000); 
get_link_list(V,1:size(V,1),{},0,'gene_network_ranklist.txt');
