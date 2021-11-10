#!/bin/bash
echo $(bc <<< "$( cat ~/hobotnica.log | wc -l) * 100 / 110")%
