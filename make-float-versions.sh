#!/bin/bash

function dtof () {
    input="$1"
    output="$2"
    gsed -E -f - "$input" > "$output" <<EOF
1i \
// Automatically generated from $(basename "$input"), DO NOT EDIT!

# Convert double to float, not anything with 'double' in the name
# that needs to be treated specially, needs to occur here
s/doubleToLongBits/floatToIntBits/g
s/double/float/g
s/Double/Float/g

# Convert back to DoubleFunction since there is no FloatFunction
s/FloatFunction/DoubleFunction/g
s/([a-z]\.apply)/(float)\1/g

# These math members are doubles, we need to cast (not abs() is missing from this list)
s/Math\.sqrt/(float)Math.sqrt/g
s/Math\.cos/(float)Math.cos/g
s/Math\.sin/(float)Math.sin/g
s/Math\.acos/(float)Math.acos/g
s/Math\.PI/(float)Math.PI/g

# Handle precision with a standard constant name
s/EPSILON = 1e-9/EPSILON = 1e-6/g

# Convert literals to have 'f' suffix
s/([0-9]+\.[0-9]+)/\1f/g
s/([0-9]+e[-+]?[0-9]+)/\1f/g
s/([a-z][0-9])d/\1f/g

# This is for hashCode implemention.  The double versions accumulate a long
# then xor the top and bottom 32 bits.  For float versions we need to remove
# the xor of the top 32 bits.
s/long h =/int h =/
s/\(int\)//g
s/\(h >>> 32\) \^ //

EOF
}

# Convert the interfaces
for file in Matrix Vector ; do
    dtof "src/main/java/edu/unc/cs/robotics/math/Double$file.java" \
         "src/main/java/edu/unc/cs/robotics/math/Float$file.java"
done

# Convert the classes and their tests
for file in AffineTransform3 Matrix3 Vec2 Vec3 Quaternion4 ; do
    main="src/main/java/edu/unc/cs/robotics/math/$file"
    test="src/test/java/edu/unc/cs/robotics/math/$file"
    dtof "${main}d.java" "${main}f.java"
    dtof "${test}dTest.java" "${test}fTest.java"
done
