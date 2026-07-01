#!/bin/bash

old_prefix="Opti_Eth_NoACT_Match_Jet_NoACT_Eth_Ret_"
new_prefix="Local_Match_Search_Att1_"

# Only change sub folders under current directory (Option 1)
# for f in "${old_prefix}"*; do
#     suffix="${f#$old_prefix}"
#     mv "$f" "${new_prefix}${suffix}"
# done

# Change current folders and sub folders (Option 2)
# find . -depth -name "${old_prefix}*" | while IFS= read -r f; do
#     dir=$(dirname "$f")
#     base=$(basename "$f")
#     suffix="${base#$old_prefix}"
#     mv "$f" "$dir/${new_prefix}${suffix}"
# done