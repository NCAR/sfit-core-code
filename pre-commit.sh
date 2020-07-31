#!/bin/bash
echo 'Version number updated(y/n)?'
while [ 1 == 1 ]; do
read text;
sleep 5;
if [[ $text == 'n' ]]; then exit 1; else exit 0; fi;
done

