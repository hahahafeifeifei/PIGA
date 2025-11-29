#!/bin/bash
# Download and build tools needed for PIGA.  They are
# - Margin
# - Secphase
# - PanGenie
# - Beagle utilities
# - MECAT2
# - Merqury
# - DeepVariant container


set -beEu -o pipefail

mainDir=$PWD
buildDir=$PWD/piga_tools_build
binDir=$CONDA_PREFIX/bin

rm -rf ${buildDir}
mkdir -p ${buildDir}
mkdir -p ${binDir}
cd ${buildDir}

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH


# Build Margin
git clone https://github.com/UCSC-nanopore-cgl/margin.git --recursive
cd margin
export LDFLAGS="-L$CONDA_PREFIX/lib -ldeflate"
cmake .
make -j 8
mv ./margin ${binDir}
mv params/phase/allParams.phase_vcf.ont.sv.json ${mainDir}/config/margin.phase_sv.json
cd ${buildDir}


# Build Secphase
git clone https://github.com/samtools/htslib.git --recursive
cd htslib
autoreconf -i
./configure --prefix=$PWD
make -j 8
make install
cd ${buildDir}

git clone https://github.com/benedictpaten/sonLib.git
cd sonLib
make -j 8
cd ${buildDir}

git clone https://github.com/mobinasri/secphase.git
cd secphase/programs/
sed -i "s/\/home\/apps\/sonLib\/lib\/sonLib.a/$(echo ${buildDir} | sed 's/\//\\\//g')\/sonLib\/lib\/sonLib.a/g" Makefile
sed -i "s/\/home\/apps\/sonLib\/C\/inc\//$(echo ${buildDir} | sed 's/\//\\\//g')\/sonLib\/C\/inc\//g" Makefile
sed -i "s/\/usr\/local\/include\/htslib/$(echo ${buildDir} | sed 's/\//\\\//g')\/htslib\/include\/htslib\//g" Makefile
make
mv bin/correct_bam ${binDir}
mv bin/secphase ${binDir}
cd ${buildDir}


# Build MECAT2
git clone https://github.com/xiaochuanle/MECAT2.git
cd MECAT2
make -j 8
mv */bin/mecat2map ${binDir}
mv */bin/mecat2pm4 ${binDir}
mv */bin/mecat2lcr ${binDir}
mv */bin/mecat2splitreads ${binDir}
mv */bin/mecat2trimbases ${binDir}
cd ${buildDir}


# Build PanGenie
git clone https://github.com/eblerjana/pangenie.git
cd pangenie
cmake .
make -j 8
mv src/PanGenie* ${binDir}
cd ${buildDir}


# Build cactus-gfa-tools
git clone https://github.com/ComparativeGenomicsToolkit/cactus-gfa-tools.git
cd cactus-gfa-tools
make -j 8
mv gaf2paf ${binDir}
mv gaffilter ${binDir}


# Download Merqury
wget https://github.com/marbl/merqury/archive/v1.3.tar.gz
tar -zxvf v1.3.tar.gz
cd merqury-1.3
sed 's/$MERQURY\/eval\///g' ~/miniforge3/envs/test4/bin/merqury.sh
mv merqury.sh ${binDir}
mv trio ${binDir}
mv util ${binDir}
mv eval ${binDir}


# Download Beagle utilities
wget https://faculty.washington.edu/browning/beagle_utilities/splitvcf.jar
chmod +x splitvcf.jar
mv ./splitvcf.jar ${binDir}

wget https://faculty.washington.edu/browning/beagle_utilities/mergevcf.jar
chmod +x mergevcf.jar
mv ./mergevcf.jar ${binDir}
cd ${buildDir}


# Download DeepVariant container
singularity pull docker://google/deepvariant
mv deepvariant*.sif ${mainDir}/config/deepvariant.sif
cd ${mainDir}
rm -rf ${buildDir}