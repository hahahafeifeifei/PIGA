#!/bin/bash
# Download and build tools needed for PIGA.  They are
# - Margin
# - Secphase
# - cactus-gfa-tools
# - orca
# - Beagle utilities
# - MECAT2
# - Merqury
# - DeepVariant

mainDir=$PWD
buildDir=$PWD/piga_tools_build
calllrDir=$PWD/scripts/call_lr_snv
inferDir=$PWD/scripts/infer_diploid_path
binDir=$CONDA_PREFIX/bin
thread=$(nproc | awk '{print ($1>8)?8:$1}')

rm -rf ${buildDir}
mkdir -p ${buildDir}
mkdir -p ${binDir}
mkdir -p ${calllrDir}
cd ${buildDir}

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
export LIBRARY_PATH=$CONDA_PREFIX/lib
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
export C_INCLUDE_PATH=$CONDA_PREFIX/include
export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include

# Build Margin
git clone https://github.com/UCSC-nanopore-cgl/margin.git --recursive
cd margin
export LDFLAGS="-L$CONDA_PREFIX/lib -ldeflate"
cmake .
make -j ${thread}
mv ./margin ${binDir}
mv params/phase/allParams.phase_vcf.ont.sv.json ${inferDir}/margin.phase_sv.json
mv params/base_params.json ${mainDir}/scripts/base_params.json
cd ${buildDir}

# Build Secphase
git clone https://github.com/samtools/htslib.git --recursive
cd htslib
autoreconf -i
./configure --prefix=$PWD
make -j ${thread}
make install
cd ${buildDir}

git clone https://github.com/benedictpaten/sonLib.git
cd sonLib
make -j ${thread}
cd ${buildDir}

git clone https://github.com/mobinasri/secphase.git
cd secphase/programs/
g++ -shared -fPIC -o submodules/edlib/libedlib.so submodules/edlib/edlib.cpp
cp submodules/edlib/libedlib.so $CONDA_PREFIX/lib
cp submodules/edlib/edlib.h $CONDA_PREFIX/include
sed -i "s/\/home\/apps\/sonLib\/lib\/sonLib.a/$(echo ${buildDir} | sed 's/\//\\\//g')\/sonLib\/lib\/sonLib.a/g" Makefile
sed -i "s/\/home\/apps\/sonLib\/C\/inc\//$(echo ${buildDir} | sed 's/\//\\\//g')\/sonLib\/C\/inc\//g" Makefile
sed -i "s/\/usr\/local\/include\/htslib/$(echo ${buildDir} | sed 's/\//\\\//g')\/htslib\/include\/htslib\//g" Makefile
sed -i "s/-Lsubmodules\/edlib/-Lsubmodules\/edlib -L$(echo ${CONDA_PREFIX}/lib | sed 's/\//\\\//g') -lstdc++/g" Makefile
make
mv bin/correct_bam ${binDir}
mv bin/secphase ${binDir}
cd ${buildDir}

# Build MECAT2
git clone https://github.com/xiaochuanle/MECAT2.git
cd MECAT2
make -j ${thread}
mv */bin/mecat2map ${binDir}
mv */bin/mecat2pm4 ${binDir}
mv */bin/mecat2lcr ${binDir}
mv */bin/mecat2splitreads ${binDir}
mv */bin/mecat2trimbases ${binDir}
cd ${buildDir}

# Build cactus-gfa-tools
git clone https://github.com/ComparativeGenomicsToolkit/cactus-gfa-tools.git
cd cactus-gfa-tools
make -j ${thread}
mv gaf2paf ${binDir}
mv gaffilter ${binDir}
cd ${buildDir}

# Build ORCA
git clone https://github.com/thocevar/orca.git
cd orca
g++ -O2 -std=c++11 -o orca.exe orca.cpp
mv orca.exe ${binDir}
cd ${buildDir}

# Download Merqury
wget https://github.com/marbl/merqury/archive/v1.3.tar.gz
tar -zxvf v1.3.tar.gz
cd merqury-1.3
mv merqury.sh ${binDir}
mv build ${binDir}
mv eval ${binDir}
mv plot ${binDir}
mv trio ${binDir}
mv util ${binDir}
cd ${buildDir}

# Download Beagle utilities
wget https://faculty.washington.edu/browning/beagle_utilities/splitvcf.jar
chmod +x splitvcf.jar
mv splitvcf.jar ${calllrDir}

wget https://faculty.washington.edu/browning/beagle_utilities/mergevcf.jar
chmod +x mergevcf.jar
mv mergevcf.jar ${calllrDir}

# Download DeepVariant container
if singularity pull docker://google/deepvariant; then
    echo "Download deepvariant successful."
else
    echo "Primary URL failed, trying backup URL."
    if singularity pull docker://docker.1ms.run/google/deepvariant; then
        echo "Download deepvariant successful."
    else
        echo "Download deepvariant failed."
    fi
fi

mv deepvariant*.sif ${calllrDir}/deepvariant.sif
cd ${mainDir}
rm -rf ${buildDir}
