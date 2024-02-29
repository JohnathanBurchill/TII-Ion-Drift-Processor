#!/bin/sh

(cd build && make) && rm -rf /tmp/0401 && ./build/tiict0401 A 2014 5 1 0302 0401 /storage/EfiCalCdfs /storage/LP /tmp
