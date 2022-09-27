FROM gitpod/workspace-full

RUN brew upgrade && brew install \
    R \
    gcc@5 \
    libtiff