#!/bin/bash
echo ""
echo "#---------- $0 (QGL) start..."
echo ""
#----- Simple script to capture basic documentation rendering logic.

echo "#-----Building QGL docs via ~:"

export CMD="mkdocs build"
echo $CMD
echo ""
$CMD

pushd .

cd site
echo ""
echo "Target QGL documents list as follows:"
pwd
ls -l


echo ""
echo "#---------- $0 (QGL) stop."
echo ""

