rm -rf build/
cmake -S bindings/src/AdjointFTR/ -B build
cmake --build build -j8