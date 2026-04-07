#!/bin/bash

echo -e "SAMN\tcollection_date" > biosample_dates.tsv

for i in $(seq 11950081 11950286); do
    SAMN="SAMN${i}"
    
    DATE=$(esearch -db biosample -query ${SAMN} | \
        efetch -format xml | \
	xtract -pattern Attribute -if "@attribute_name" -equals "collection_date" -element .)

    echo "${SAMN} = ${DATE}"
    echo -e "${SAMN}\t${DATE}" >> biosample_dates.tsv
done

sed 's/Attribute "\(.*\)"/\1/' biosample_dates.tsv > cleaned_dates.tsv