FROM gitpod/workspace-full

RUN brew upgrade && brew install \
    R \
    gcc@5 \
    libtiff

FROM gitpod/workspace-full
RUN pip3 install radian