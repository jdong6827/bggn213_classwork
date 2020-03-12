Structure vased drug design
================

## Download and process starting strucure

Here we download and clean up the HIV-Pr strucure(PDB code: 1HSG) from
the main PDB database. We will make a seperate set of “protein-only” and
“ligand only” PDB files.

``` r
library("bio3d")
file.name<-get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

We will use `read.pdb()`, `atom.select()` and `write.pdb()` functions to
make our seperate “protein-only” and “ligand-only” PDB files.

``` r
pdb<-read.pdb(file.name)
```

``` r
prot<-atom.select(pdb,"protein", value=TRUE)
ligd<-atom.select(pdb,"ligand", value=TRUE)
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(ligd, file="1hsg_ligand.pdb")
```

## Read docking results

Read in the output of docking and make a PDB file for viewing in VMD or
PyMol.

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```
