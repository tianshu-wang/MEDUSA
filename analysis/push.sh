#!/bin/bash

my_dir="$(dirname "$0")"
source $my_dir/dirlist.sh

for d in ${dirlist[@]}; do
  mkdir -p -- ../$d/analysis
  echo "Pushing" $1 "to" $d "..."
  cp $1 ../$d/analysis
done
