#!/bin/bash

old_prefix="Opti_Eth_NoACT_Match_Jet_NoACT_Eth_Ret_"
new_prefix="Local_Match_Search_Att1_"

for f in "${old_prefix}"*; do
    suffix="${f#$old_prefix}"
    mv "$f" "${new_prefix}${suffix}"
done