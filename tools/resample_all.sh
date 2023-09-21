#!/bin/bash

path="/home/jradema/Documents/sk/oceans/rx_outputs/River_PAB2_Van_Atta_Range_0124_0127_0130"

ls -d /home/jradema/Documents/sk/oceans/rx_outputs/VAB_Elements/* | xargs -I{} -P15 python3 convert_sig.py {}