FROM gitpod/workspace-full

RUN brew upgrade && brew install -y \
    R \
    gcc@5 \
    libtiff