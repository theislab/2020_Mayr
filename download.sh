#!/bin/bash

# downloads 9Gb dataset files
wget https://drive.google.com/file/d/13vf6Fcy6cCJUuGvbnj5sQDhayLRq7op1 -O mayr_et_al.tar.gz
mkdir data; cd data
gunzip ../mayr_et_al.tar.gz
tar -xvf ../mayr_et_al.tar
