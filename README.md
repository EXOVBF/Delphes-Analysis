Delphes-Analysis
================

 * downloading and compiling Delphes on the SLC6 from lxplus
```
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.0.12.tar.gz
tar -zxf Delphes-3.0.12.tar.gz
cd Delphes-3.0.12/
cmsrel CMSSW_6_2_0
cd CMSSW_6_2_0
cmsenv
cd ../
source ../Decay/setup_slc6.sh
make -j 16
```
