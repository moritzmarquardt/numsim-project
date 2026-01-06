#!/bin/bash

# Configuration
TEMPLATE="lid_driven_cavity.txt.template"
EXECUTABLE="./build/numsim"

# Use seq to generate 101 values
# Re: 500 down to 150 (decrementing by 3.5 per step approx)
# Velocity: 0.5 to 1.5 (incrementing by 0.01 per step)
RE_VALUES=$(seq 500 10 1500 | head -n 101)
VEL_VALUES=$(seq 0.5 0.01 1.5 | head -n 101)

# Convert strings to arrays to iterate easily
re_arr=($RE_VALUES)
vel_arr=($VEL_VALUES)

# Loop through the indices
for i in {0..100}; do
    RE=${re_arr[$i]}
    VEL=${vel_arr[$i]}
    
    # Define a unique name for this run's config file
    CURRENT_CONFIG="config_re${RE}_v${VEL}.txt"
    
    echo "Running simulation $i: Re=$RE, Velocity=$VEL"
    
    # Use sed to replace placeholders and create a temporary config file
    sed -e "s/RE_VALUE/$RE/" \
        -e "s/VELOCITY_VALUE/$VEL/" "$TEMPLATE" > "$CURRENT_CONFIG"
    
    # Run the executable with the generated config
    $EXECUTABLE "$CURRENT_CONFIG"
    
    # Rename output to include simulation index to avoid overwrite
    mv ./out/output_0000.vti ./out/output_${i}.vti

    # Optional: remove the temp config file after the run to keep directory clean
    rm "$CURRENT_CONFIG"
done