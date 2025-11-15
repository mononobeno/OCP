FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

# 基本ツール & MPI & ビルド環境
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential cmake git wget ca-certificates \
      openmpi-bin libopenmpi-dev \
      libfftw3-dev && \
    rm -rf /var/lib/apt/lists/*

# GROMACS 2023.x をソースからビルド (CUDA + MPI)
WORKDIR /opt
ENV GMX_VERSION=2023.3

RUN wget https://ftp.gromacs.org/gromacs/gromacs-${GMX_VERSION}.tar.gz && \
    tar xfz gromacs-${GMX_VERSION}.tar.gz && \
    rm gromacs-${GMX_VERSION}.tar.gz && \
    mkdir -p gromacs-${GMX_VERSION}-build && \
    cd gromacs-${GMX_VERSION}-build && \
    cmake ../gromacs-${GMX_VERSION} \
      -DGMX_BUILD_OWN_FFTW=ON \
      -DGMX_MPI=ON \
      -DGMX_GPU=CUDA \
      -DGMX_USE_RDTSCP=OFF \
      -DREGRESSIONTEST_DOWNLOAD=OFF \
      -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs && \
    make -j$(nproc) && \
    make install && \
    rm -rf /opt/gromacs-${GMX_VERSION} /opt/gromacs-${GMX_VERSION}-build

# 環境変数設定
ENV PATH=/usr/local/gromacs/bin:$PATH
ENV GMXBIN=/usr/local/gromacs/bin
ENV GMXLIB=/usr/local/gromacs/share/gromacs
ENV GMX_MAXBACKUP=-1

# デフォルトは gmx_mpi のバージョンを表示
CMD ["gmx_mpi", "--version"]
