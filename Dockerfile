FROM ubuntu:18.04 as builder
RUN apt update && apt upgrade -y && \
    apt install -y cmake g++ make nano libeigen3-dev
RUN cp -r /usr/include/eigen3/Eigen /usr/include
ADD . /opt/sources
WORKDIR /opt/sources
RUN mkdir build && cd build && cmake .. && make && make install

FROM ubuntu:18.04
WORKDIR /usr/bin
COPY --from=builder /usr/local/bin .
RUN cp ./opendlv-sim-can-rhino ./opendlv-sim-can-rhino-2
