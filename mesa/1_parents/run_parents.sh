#!/bin/bash
set -e
for m in 0.7 0.6; do
  cp inlist_parent_$m inlist_project
  ./rn
done
