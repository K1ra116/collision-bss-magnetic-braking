#!/bin/bash
set -e
cp inlist_relax inlist_project
./rn
if [ "$1" == "evolve" ]; then
  cp inlist_evolve inlist_project
  ./rn
fi
