# math

This is a collection of mathematical classes useful in my robotics projects.  It contains 3D transforms, rotations, vectors, and some miscellaneous math methods.

### Matrix and Vector Usage Patterns
There are API patterns to the matrix and vector classes in library: results pre-allocation and method chaining.

#### Pre-allocation of Results
Pre-allocate results, and use the result as the object of the method call.   The pre-allocation allows the caller to completely control object allocation and thus reuse objects.  This can (usually) avoid excessive allocation and garbage collection as a result of mathematical routines.  Also by having the result object be the base of the call, it follows mathematical conventions of writing "x = a * b", for example:

```java
Vec3d x = new Vec3d(); // preallocate result
x.cross(a, b); // compute cross product of a and b and store result in x
```
#### Method Chaining

Where possible, methods return `this`.  This allows for method chaining, where the result from the first operation is used immediately as an operator in the second.  For example:
```java
Vec3d v = new Vec3d();
v.sub(a,b).mul(0.5).add(c); // compute v = (a-b) * 0.5 + c
```

### Double / Float Versions

This code is primarily written in `double` precision math.  The `make-float-versions.sh` script converts the double versions to float version.  This is a very simple sed-script based conversion which relies upon a few rules:

* For all floating-point values, use a decimal, e.g., instead of "0", use "0.0".  The script will append "f" to make them "0.0f".
* Follow the existing classes for computing hash codes based upon doubles, in particular, make sure to accumulate a long, and then return the XOR'd value in the same form.


License
----

BSD 2-Clause License
