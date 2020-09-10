#!/bin/bash

echo "start"

R CMD BATCH tcar1.R
echo "tcar1 finished"

R CMD BATCH tcar2.R
echo "tcar2 finished"

R CMD BATCH tcarw.R
echo "tcarw finished"

R CMD BATCH tcare1.R
echo "tcare1 finished"

R CMD BATCH tcare2.R
echo "tcare2 finished"
