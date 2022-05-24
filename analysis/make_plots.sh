#!/bin/bash

cwd=`pwd`
my_dir="$(dirname "$0")"
source $my_dir/dirlist.sh

for d in ${dirlist[@]}
do
  cd ../$d/analysis
  python $*
  cd $cwd
done

