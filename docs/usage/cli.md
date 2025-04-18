TIES can be access via both command line and python interface.

In the smallest example one can carry out a superimposition 
employing only two ligands
```shell
ties --ligands l03.mol2 l02.mol2
```

Ideally these .mol2 files already have a 
pre-assigned charges in the last column for each atom. 
See for example [MCL1 case](https://github.com/UCL-CCS/TIES/blob/master/ties/examples/mol2_2ligands_MCL1/l02.mol2).

In the case of this example, we are working on molecules 
that are negatively charges (-1e), and we need to specify it:
```shell
ties --ligands l03.mol2 l02.mol2 --ligand-net-charge -1
```

The order the of the ligands matters and more ligands can be passed. 
This command creates by default a `ties-input` directory with all output 
files. These include `meta_*_*.json` files which contain the 
details about how the ligands were superimposed, and what 
configuration was used. The general directory structure 
will look like this:

```shell
    ties
    ├── mol2
    │    ├── l02
    │    └── l03
    ├── prep
    │   ├── ligand_frcmods
    │   │   ├── l02
    │   │   └── l03
    │   └── morph_frcmods
    │       └── tests
    │           └── l02_l03
    └── ties-l02-l03
        └── lig
            └── build
```

Note that all the output generated by ambertools is stored, 
and can be investigated.

The full RBFE requires also the protein, as well 
as the net charge of the ligands used in the transformation:
```shell
ties -l l02.mol2 l03.mol2 -nc -1 --protein protein.pdb
```

Check all the options with
```shell
ties -h
```


!!! warning 
    This code is currently experimental and under active development.
    If you notice any problems, report them. *