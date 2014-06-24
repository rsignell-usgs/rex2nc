#!/bin/sh

# Convert from 3 letter slosh domains to 4 letter domains or vice versa
# If the 3 letter domain is a "p" domain, the p is not output.

INPUT=${1}
if [ x"${INPUT}" = x"" ]; then
    echo "Input basin (3 or 4 characters) must be provided on the command line."
    exit -1
fi

BNT_FILE="sloshdsp.bnt"

# See how many characters we are starting with
NUM_CHAR=`echo ${INPUT} | wc -m`
NUM_CHAR=`expr ${NUM_CHAR} - 1 `

# Set the field seperator to :
OFS="${IFS}"
IFS=':'

# If 3 characters then convert to 4
if [ x"${NUM_CHAR}" = x"3" ]; then
    while read LINE; do
	
       ARRAY_LINE=($LINE)  # Uses IFS to seperate
       
       FIELD2=${ARRAY_LINE[1]}
       FIELD3=${ARRAY_LINE[2]}
       FIELD4=${ARRAY_LINE[3]}

       if [ x"${FIELD4}" = x"${INPUT}" ]; then
	   if [ x"${FIELD2}" = x"p" ]; then
#              Drop the "p"
	       echo "${FIELD3}"
	   else
	       echo "${FIELD2}${FIELD3}"
	   fi
       fi

    done < ${BNT_FILE}
fi

# If 4 characters then convert to 3
if [ x"${NUM_CHAR}" = x"4" ]; then

    while read LINE; do

       ARRAY_LINE=($LINE)  # Uses IFS to seperate
       
       FIELD2=${ARRAY_LINE[1]}
       FIELD3=${ARRAY_LINE[2]}
       FIELD4=${ARRAY_LINE[3]}

       if [ x"${FIELD2}${FIELD3}" = x"${INPUT}" ]; then
	   echo "${FIELD4}"
       fi

    done < ${BNT_FILE}

fi

# Restore the field seperator
IFS="${OFS}"
