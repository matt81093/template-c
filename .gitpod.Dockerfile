FROM gitpod/workspace-full

RUN sudo apt-get update && \
    sudo apt-get install -y \ 
    valgrind libcmocka-dev cmocka-doc libcmocka0 \
    openmpi-bin openmpi-dev openmpi-common openmpi-doc libopenmpi-dev

