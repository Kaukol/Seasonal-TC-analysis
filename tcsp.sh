#!/bin/bash

echo "start"

R CMD BATCH tcsp.R
echo "tcar1 finished"

R CMD BATCH tcspe1.R
echo "tcar2 finished"

R CMD BATCH tcspe2.R
echo "tcarw finished"

R CMD BATCH tcspw1.R
echo "tcare1 finished"

R CMD BATCH tcspw2.R
echo "tcare2 finished"
