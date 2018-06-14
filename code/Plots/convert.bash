#!/bin/bash          

for i in *.eps; do
    [ -f "$i" ] || break
    epstopdf $i
done

cp *.pdf ~/Dropbox/=AOAS/revision/plots/
cp *.pdf ~/Dropbox/=AOAS/latex/plots/
