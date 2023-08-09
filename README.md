gloptMP stands for global optimization multi-process

Compilation:
```
mkdir gloptMP/build
cd gloptMP/build
cmake ..
make
```

Usage examples:
```
./gloptMP --help
```
Output help message

```
./gloptMP -d 2 -v 50 -b 8 -i 1000 -p 4 -e 1e-3
```
Optimize functional for 2D, starting with 50 points, with 8 bit resolution for each dimension, perform 1000 iterations, adding 4 points at once. Stopping criterion is 1e-3 (a small one)
