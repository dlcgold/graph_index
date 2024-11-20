```bash
git clone --recursive https://github.com/dlcgold/graph_index.git
cd graph_index/gfatools
make -j2
cd ../mfmi
make -j2
cd ..
make -j2

./gindex index -i example/test.gfa -o example/test-index -t 5 - c 5
./gindex query -i example/test-index -q example/test.fq
# Matches (as nodes in .gfa):
# @read1    6
# @read2    4
# @read2    4
```
