#!/bin/sh

(cd build && make) && rm -rf /tmp/0402 && ./build/tiict0402 A 2014 5 1 0302 0402 /storage/EfiCalCdfs /storage/LP /tmp
