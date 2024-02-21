* [Requirements](#requirements)
* [Macros](#macros)
  * [SelectClusters.C](#selectclustersc)
  * [ClusterCOG.C](#clustercogc)
  * [CompareClusters.C](#compareclustersc)
* [Utils](#utils)
  * [CCDBUtils.h](#ccdbutilsh)
  * [DataUtils.h](#datautilsh)
  * [PreClusterUtils.h](#preclusterutilsh)
  * [ClusterUtils.h](#clusterutilsh)

# Requirements

- The reconstruction (clustering) must have been run with the option `--attach-initial-precluster` so that all the digits from the original precluster are associated to the reconstructed cluster(s).
- The reconstruction (tracking) must have been run without the option `--digits` so that the indices of the digits associated to the clusters associated to the tracks in mchtracks.root point to the digits stored in mchclusters.root.

These 2 requirements are needed to be able to retrieve the number of clusters reconstructed from a given precluster and select isolated clusters.

# Macros

## SelectClusters.C

### Description
- Loop over muon (MCH+MID) tracks and refit them to get the track parameters at each cluster.
- If `applyTrackSelection = true`, extrapolate the tracks to the vertex and apply the standard track selections (eta, Rabs, pDCA) + p > 10 GeV/c.
- Select the isolated clusters (i.e. cases where the clustering found only one cluster in the precluster), with digits on both cathodes.

### Output
The ouput is a root file (default = `clusters.root`) containing:
* a tree named `data`. Each entry contains:
    * a selected cluster: `o2::mch::Cluster`
    * the associated track parameters: `o2::mch::TrackParamStruct`
    * the associated precluster (= a vector of associated digits): `std::vector<o2::mch::Digit>`
* canvases with some control plots showing the characteristics of the selected preclusters

### Example
```shell
gSystem->Load("libO2MCHMappingImpl4")
.x SelectClusters.C+(529691,"mchclusters.root","mchtracks.root","muontracks.root",true)
```

## ClusterCOG.C

### Description
- Loop over the clusters selected by the macro [SelectClusters.C](#selectclustersc).
- Apply further selections on e.g. the size of the associated preclusters (= number of associated digits) or their charge asymmetry between the two cathodes.
- Create new clusters positioned at the center-of-gravity (COG) of the associated digits.

### Output
The ouput is a root file (default = `newclusters.root`) containing:
* a tree named `data`. Each entry contains:
    * a selected cluster: `o2::mch::Cluster`
    * the associated track parameters: `o2::mch::TrackParamStruct`
    * the associated precluster (= a vector of associated digits): `std::vector<o2::mch::Digit>`
    * the corresponding new cluster positioned at the COG: `o2::mch::Cluster`
* canvases with some control plots showing the characteristics of the selected preclusters

### Example
```shell
gSystem->Load("libO2MCHMappingImpl4")
.x ClusterCOG.C+(529691)
```

## CompareClusters.C

### Function 1
Compare the original and the new clusters stored in the same file (produced by e.g. [ClusterCOG.C](#clustercogc)):
- Draw the residuals between the original and new cluster positions
- Draw the residuals between the original(new) cluster and the track positions
- Fit the cluster - track residuals with a Crystal-Ball to extract a pseudo-resolution
- Draw the ratio of the cluster - track residuals between the original and the new clusters

```shell
.x CompareClusters.C+
```

### Function 2
Compare either the original clusters or the new clusters stored in two different files (produced by e.g. [SelectClusters.C](#selectclustersc) or [ClusterCOG.C](#clustercogc)):
- Draw the cluster - track residuals from the two files
- Fit the cluster - track residuals with a Crystal-Ball to extract a pseudo-resolution
- Draw the ratio of the cluster - track residuals between the two files

```shell
.x CompareClusters.C+("newclusters1.root","newclusters2.root",true)
```

# Utils

## CCDBUtils.h

Utility function to load geometry and/or magnetic field from CCDB or local snapshot, and prepare track extrapolation if needed.

## DataUtils.h

Utility function to open a root file and get a tree.

## PreClusterUtils.h

- Contains utility functions to determine the characteristics of a precluster:
  - mono-cathode (all digits on the same cathode) or not
  - size (= number of digits) on each cathode
  - total charge on each cathode
  - center-of-gravity (2 methods)
- Contains utility functions to produce, fill and draw some control plots.

To be used, they require the MCH mapping to be loaded (`gSystem->Load("libO2MCHMappingImpl4")`).

## ClusterUtils.h

Utility function to produce a new cluster in the global coordinate system from a given position in the local coordinate system of a DE. The DE ID is extracted from the provided cluster uid.

To be used, it requires the geometry to be loaded (e.g. from the CCDB).
