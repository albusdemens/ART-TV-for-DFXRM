# ART-TV

Code to reconstruct the shape of a grain from a topo-tomography dataset collected using dark-field X-ray microscopy (DFXRM) at beamline BL06 of the European Synchrotron Radiation Facility (ESRF).

The reconstruction algorithm combines the Algebraic Reconstruction Technique (ART) and Total Variation (TV) <sup>[1](#myfootnote1)</sup>. In the reconstruction, orientation of the voxels is calculated from the acquired frames using [recon3d](https://github.com/albusdemens/Recon3D).

Contributors: A.Cereser, M. Busi

<a name="myfootnote1">1</a>: LaRoque, S. J., Sidky, E. Y., & Pan, X. (2008). Accurate image reconstruction from few-view and limited-angle data in diffraction tomography. JOSA A, 25(7), 1772-1782.
