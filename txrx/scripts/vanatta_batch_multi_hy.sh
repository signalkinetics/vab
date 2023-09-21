#!/bin/bash

dists="90 110 130 150"


echo -e $dists | xargs -d ' ' -I{} -P2 python3 vanatta_batch_eq_multi_hyd.py {} --volt 40
