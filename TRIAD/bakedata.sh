#!/usr/bin/env bash

# bash isn't the right tool for this, maybe make a python script later

DATA_FILE=$'src/Geomag/Data.h'
DATA_HEADER=$'#ifndef _geomag_data_h\n#define _geomag_data_h\n\n#include <Geomag/GeomagLib.h>\n'
DATA_FOOTER=$'#endif'

WMM_STRUCT=$'typedef struct wmm_element_s\n{\n    double coeffg;\n    double coeffh;\n    double slopeg;\n    double slopeh;\n} wmm_element_t;\n'
WMM_DATA_HEADER=$'static const wmm_element_t geomag_wmm_elements[] =\n{'
WMM_FILE=$'WMM.COF'
WMM_DATA_FOOTER=$'};\n'

WMM_TEST_STRUCT=$'typedef struct wmm_test_element_s\n{\n    double time;\n    double alt, lat, lon;\n    double x, y, z;\n} wmm_test_element_t;\n'
WMM_TEST_DATA_HEADER=$'static const wmm_test_element_t geomag_wmm_test_elements[] =\n{'
WMM_TEST_FILE=$'WMM2025_TestValues.txt'
WMM_TEST_DATA_FOOTER=$'};\n'

rm "$DATA_FILE"

echo "$DATA_HEADER" >> "$DATA_FILE"

echo "$WMM_STRUCT" >> "$DATA_FILE"
echo "$WMM_TEST_STRUCT" >> "$DATA_FILE"

echo "$WMM_DATA_HEADER" >> "$DATA_FILE"
WMM_TERMINATOR=$'999999999999999999999999999999999999999999999999'
tail -n +2 "$WMM_FILE" | while read -r line; do
    if [ "$line" == "$WMM_TERMINATOR" ]; then
        break
    fi

    echo -n "    { " >> "$DATA_FILE"

    i=0
    for coeff in $line; do
        if [ "$i" -gt "1" ]; then
            echo -n "$coeff, " >> "$DATA_FILE"
        fi
        
        ((i++))
    done

    echo "}," >> "$DATA_FILE"
done
echo "$WMM_DATA_FOOTER" >> "$DATA_FILE"

echo "$WMM_TEST_DATA_HEADER" >> "$DATA_FILE"
tail -n +19 "$WMM_TEST_FILE" | while read -r line; do
    
    echo -n "    { " >> "$DATA_FILE"
    
    i=0
    for val in $line; do
        if [ "$i" -eq "0" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        elif [ "$i" -eq "1" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        elif [ "$i" -eq "2" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        elif [ "$i" -eq "3" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        elif [ "$i" -eq "7" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        elif [ "$i" -eq "8" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        elif [ "$i" -eq "9" ]; then
            echo -n "$val, " >> "$DATA_FILE"
        fi
        
        ((i++))
    done

    echo "}," >> "$DATA_FILE"
done
echo "$WMM_TEST_DATA_FOOTER" >> "$DATA_FILE"

echo "$DATA_FOOTER" >> "$DATA_FILE"