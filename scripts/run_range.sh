#!/bin/bash
#=============================================================================
# Production driver: process a range of years/months in sequence.
#
# Usage:
#   ./run_range.sh <YEAR_START> <YEAR_END> [MONTH_START] [MONTH_END]
#
# Examples:
#   ./run_range.sh 2014 2014               # all 12 months of 2014
#   ./run_range.sh 2014 2014 1 3           # Jan-Mar 2014
#   ./run_range.sh 2014 2016               # all 12 months of 2014, 2015, 2016
#=============================================================================
set -euo pipefail

if [ $# -lt 2 ] || [ $# -gt 4 ]; then
    echo "Usage: $0 <YEAR_START> <YEAR_END> [MONTH_START] [MONTH_END]"
    echo "       (months default to 1..12)"
    exit 1
fi

YEAR_START=$1
YEAR_END=$2
MONTH_START=${3:-1}
MONTH_END=${4:-12}

if [ "$YEAR_START" -gt "$YEAR_END" ]; then
    echo "ERROR: YEAR_START ($YEAR_START) > YEAR_END ($YEAR_END)"
    exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

TOTAL=0
DONE=0
FAILED=()
for YEAR in $(seq "$YEAR_START" "$YEAR_END"); do
    for MONTH in $(seq "$MONTH_START" "$MONTH_END"); do
        TOTAL=$((TOTAL + 1))
        echo ""
        echo "########################################################################"
        echo "  [$((DONE + 1))/$TOTAL] Processing y=$YEAR m=$MONTH"
        echo "########################################################################"
        if "$SCRIPT_DIR/run_month.sh" "$YEAR" "$MONTH"; then
            DONE=$((DONE + 1))
        else
            FAILED+=("y=$YEAR m=$MONTH")
            echo "  ✗ FAILED y=$YEAR m=$MONTH (continuing...)"
        fi
    done
done

echo ""
echo "########################################################################"
echo "  COMPLETE: $DONE/$TOTAL months processed successfully"
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "  Failures:"
    for f in "${FAILED[@]}"; do
        echo "    - $f"
    done
    exit 1
fi
echo "########################################################################"
