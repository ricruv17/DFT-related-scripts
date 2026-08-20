#!/bin/bash

# sum_xsf.sh — sum the volumetric data of two XSF files
# Usage: ./sum_xsf.sh file1.xsf file2.xsf
# Output: sum.xsf

set -e

f1="$1"
f2="$2"

# ---------------------------------------------------------------------------
# Helper: extract the line number of a pattern (first match)
# ---------------------------------------------------------------------------
grep_line() { grep -n "$1" "$2" | head -1 | cut -d: -f1; }

# ---------------------------------------------------------------------------
# Locate the BEGIN_DATAGRID line in each file — everything above it
# (CRYSTAL, PRIMVEC, PRIMCOORD, atom list) is the "prefix" we keep from f1.
# ---------------------------------------------------------------------------
datagrid_line1=$(grep_line "BEGIN_DATAGRID_3D" "$f1" | head -1)
datagrid_line2=$(grep_line "BEGIN_DATAGRID_3D" "$f2" | head -1)

# The header = everything up to and including BEGIN_DATAGRID_3D_xxx (1 line after BEGIN_BLOCK)
# We keep the full structural header from file1 unchanged.
header_end=$((datagrid_line1 + 5))   # label + nx ny nz + origin + 3 cell vectors = 5 lines after BEGIN

# ---------------------------------------------------------------------------
# Split file1 into: prefix (up to and including grid vectors) | volumetric data | suffix
# ---------------------------------------------------------------------------
prefix_lines=$header_end

vol_start1=$((prefix_lines + 1))
# Find END_DATAGRID line to know where volumetric data stops
end_line1=$(grep_line "END_DATAGRID_3D" "$f1")
vol_end1=$((end_line1 - 1))

vol_start2=$((datagrid_line2 + 6))   # same offset for file2
end_line2=$(grep_line "END_DATAGRID_3D" "$f2")
vol_end2=$((end_line2 - 1))

# Extract parts
head -n "$prefix_lines" "$f1"        > xsf_prefix
sed -n "${vol_start1},${vol_end1}p"  "$f1" > vol1
sed -n "${vol_start2},${vol_end2}p"  "$f2" > vol2
sed -n "${end_line1},\$p"            "$f1" > xsf_suffix   # END_DATAGRID onward

# ---------------------------------------------------------------------------
# Sum volumetric data element-by-element with awk
# (same paste | awk approach as the cube script)
# ---------------------------------------------------------------------------
cat > tmp_xsf.awk <<'EOF'
{
    for (i = 1; i <= NF/2; i++) {
        j = i + NF/2
        printf("%12.6e ", $i + $j)
    }
    printf("\n")
}
EOF

paste vol1 vol2 | awk -f tmp_xsf.awk > vol_sum

# ---------------------------------------------------------------------------
# Reassemble
# ---------------------------------------------------------------------------
cat xsf_prefix vol_sum xsf_suffix > sum.xsf

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
rm tmp_xsf.awk xsf_prefix vol1 vol2 xsf_suffix vol_sum

echo "Done → sum.xsf"
