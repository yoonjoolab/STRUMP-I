
FROM ubuntu:focal

RUN apt-get update && apt-get upgrade -y && \
apt install -y apt-utils build-essential python3.8-dev python3-distutils wget

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1

WORKDIR /opt
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py
RUN python3 -m pip install numpy pandas psutil freesasa




