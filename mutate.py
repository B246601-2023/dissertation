import pyvolve
mytree = pyvolve.read_tree(file = "tree_constant.tre")

#custom_mu = {"AC":0.5, "AG":0.25, "AT":1.23, "CG":0.55, "CT":1.22, "GT":0.47}
#freqs = [0.1, 0.45, 0.3, 0.15]

# define a nucleotide model 
mymodel = pyvolve.Model(
    "nucleotide",
    {"kappa":3.5},
    alpha = 0.5, num_categories = 3, pinv = 0.25
    )

my_partition  =  pyvolve.Partition ( models  =  mymodel ,  size  =  10000 )
evolver = pyvolve.Evolver(tree = mytree,partitions = my_partition)
evolver(ratefile = "custom_ratefile.txt", infofile = "custom_infofile.txt", seqfile
= "custom_seqfile.fasta" )
