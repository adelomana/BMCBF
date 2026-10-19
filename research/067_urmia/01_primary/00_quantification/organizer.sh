#!/usr/bin/env bash


for f in *_1.fq.gz; do
    label="${f%_1.fq.gz}"
    mkdir -p "$label"
    mv "${label}_1.fq.gz" "${label}_2.fq.gz" "$label/"
done

