#!/bin/sh
# Compile and run. Extra arguments: test | demo | scale
mkdir -p out
javac -d out -sourcepath src src/texthack/app/Main.java && java -Xmx2g -cp out texthack.app.Main "$@"
