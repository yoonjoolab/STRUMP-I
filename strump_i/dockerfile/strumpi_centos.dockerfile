
FROM centos:centos8

RUN yum update -y

RUN yum -y groupinstall 'Development Tools'
RUN yum -y install python38-devel wget
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1
WORKDIR /opt
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py
RUN python3 -m pip install numpy pandas psutil freesasa



