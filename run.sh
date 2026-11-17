#!/bin/sh

(cd build && make) && rm -rf /tmp/0501 && ./build/tiict0501 A 2015 10 18 0501 0501 /swarm_synology/cefi_db/EfiCalCdfs /storage/LP /tmp
