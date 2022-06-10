FROM centos:7
LABEL Author="Adrian Flannery" \
      Maintainer="aflanry@gmail.com"

WORKDIR /data/

ADD bcl-convert-3.10.5-2.el7.x86_64.rpm ./
RUN yum -y localinstall bcl-convert-3.10.5-2.el7.x86_64.rpm
RUN rm bcl-convert-3.10.5-2.el7.x86_64.rpm

# BUILD_CMD: docker build -t aflanry/bcl-convert:3.10.5 ./
