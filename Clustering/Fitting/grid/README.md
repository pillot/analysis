* [Principle](#principle)
* [Procedure](#procedure)
* [Intermediate merging](#intermediate-merging)
* [Final merging](#final-merging)
* [Local testing](#local-testing)

# Principle

Merge the input files on the grid then run the SelectClusters.C macro on them. This requires to be in an O2 environment with a valid alien token.

# Procedure

### 1. produce an xml collection of mchtracks.root files to be processed

```shell
alien_find -x - /alice/data/2024/LHC24ac/549586/apass1_muon_matching mchtracks.root > data.xml
```

This collection will be used to list the directories on the grid where the input data (mchclusters.root, mchtracks.root and muontracks.root) will be found.

### 2. copy the scripts, macros and collection to the working directory on the grid

The path to the working directory on the grid is either absolute or relative to the home directory.

```shell
path-to-CopyFiles-script/CopyFiles.sh alien://run3/data/2024/LHC24ac/549586/apass1_muon_matching/clusters_pT1GeV-p10GeV

alien_cp file:data.xml alien://run3/data/2024/LHC24ac/549586/apass1_muon_matching/clusters_pT1GeV-p10GeV/data.xml
```

The script `CopyFiles.sh` also produces and copies the .jdl files. The software package used on the grid (e.g. `VO_ALICE@O2PDPSuite::daily-20250214-0000-1`) can be optionally changed by giving it as a second parameter of the script.

### 3. submit the jobs

From the working directory on the grid:
```shell
submit SelectClusters.jdl 549586
```
The run number is necessary for SelectCluster.C to setup the geometry and magnetic field.

The outputs will be store in the subdirectory `results/` of the working directory on the grid.

# Intermediate merging

For large number of output files, intermediate merging steps can be run on the grid. Similar procedure as above must be followed for each merging stage `i`:
1. produce an xml collection `stagei.xml` of outputs to be merged
2. copy it to the working directory on the grid
3. submit the merging jobs: `submit merge.jdl i`

The merged outputs will be store in the subdirectory `stagei/` of the working directory on the grid.

# Final merging

The final merging in the local working directory is done following these steps:
1. produce an xml collection `results.xml` of outputs to be merged
2. copy the macro `merge.C` to the local working directory
3. run the merging script locally: `path-to-merge-script/merge.sh results.xml`

# Local testing

It is possible to test the code locally with an xml collection (produced as described in the step 1. of the [procedure](#procedure), limiting the number of inputs with the option `-l <count>`), or with a text file containing a list of local dir/file.

The scripts and macros can be copied to the test directory with:
```shell
path-to-CopyFiles-script/CopyFiles.sh testdir
```

Then the code can be run locally from the test directory with:
```shell
./SelectClusters.sh data.xml 549586
```
