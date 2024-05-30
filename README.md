```bash
git clone --recursive https://github.com/dlcgold/graph_index.git
cd graph_index/gfatools
make -j2
cd ../mfmi
make -j2
cd ..
make -j2

./gindex -i example/test.gfa -o example/test-index
./gindexquery -i example/test-index -q example/test.fq
# Matches (as nodes in .gfa):
# 1: 6, 6
# 2: 1, 5, 6, 6, 7 
# 3: 1, 2, 3
# 4: None
```