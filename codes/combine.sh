#!/bin/bash

# set the path
folder1=../0partitions
folder2=../1stability

convert +append $folder1/$1_subunit_*.png $folder2/$1.png $1.jpeg

