#!/bin/sh

(cd build && make) && rm -rf /tmp/0501 && ./build/tiict0501 A 2014 5 1 0302 0501 /storage/EfiCalCdfs /storage/LP /tmp
