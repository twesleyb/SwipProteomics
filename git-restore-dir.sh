#!/usr/bin/env bash

# restore a directory from a previous commit
HASH="850b2a766b8271d377abfda893ded27c58792599"

git checkout "$HASH" -- ~/SwipProteomics/refs/

# restore a directory from the previous commit
#git checkout HEAD~1 -- path/to/the/folder/ 
