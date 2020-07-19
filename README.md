# gl-math

A small math library aimed at gamedev that provides 4x4 float matrix, vector, and quaternion operations. This fork implements the library as a standard R6RS scheme library.

## Documentation
gl-math provides a number of functions for working with 4x4 matrices (plus a handful of others). The functionality is similar to what can be found in the [glm egg](http://wiki.call-cc.org/eggref/4/glm).

gl-math expects matrices, vectors, and quaternions to be vectors. The vector must be 16 elements long, 3 elements long, or 4 elements long for matrices, vectors, or quaternions, respectively.

gl-math operates on matrices in a column-major fashion in correspondence with OpenGL (e.g. translation components are at indices 12, 13, and 14). Vectors are arranged as (`(x y z)`), and quaternions as (`(x y z w)`).

### Matrix operations
    [procedure] (m* A B)

Multiply matrix `A` with matrix `B`.

    [procedure] (m*s A S)

Multiply matrix `A` by scalar `S`.

    [procedure] (m+ A B)

Add matrix `A` by matrix `B`.

    [procedure] (m- A B)

Subtract matrix `B` from matrix `A`.

    [procedure] (mat4-identity)

Return an identity matrix.

    [procedure] (translation VECTOR)

Return the translation matrix given by `VECTOR`.

    [procedure] (translate VECTOR MATRIX)

Translate `MATRIX` by `VECTOR`.

    [procedure] (x-rotation ANGLE)

Return the rotation matrix given by a rotation of `ANGLE` radians around the x-axis.

    [procedure] (rotate-x ANGLE MATRIX)

Rotate `MATRIX` around the x-axis by `ANGLE` radians.

    [procedure] (y-rotation ANGLE)

Return the rotation matrix given by a rotation of `ANGLE` radians around the y-axis.

    [procedure] (rotate-y ANGLE MATRIX)

Rotate `MATRIX` around the y-axis by `ANGLE` radians.

    [procedure] (z-rotation ANGLE)

Return the rotation matrix given by a rotation of `ANGLE` radians around the z-axis.

    [procedure] (rotate-z ANGLE MATRIX)

Rotate `MATRIX` around the z-axis by `ANGLE` radians.

    [procedure] (axis-angle-rotation AXIS ANGLE)

Return the rotation matrix given by a rotation of `ANGLE` radians around the vector `AXIS`.

    [procedure] (rotate-axis-angle AXIS ANGLE MATRIX)

Rotate `MATRIX` around the vector `AXIS` by `ANGLE` radians.

    [procedure] (quaternion-rotation Q) 

Return the rotation matrix given by the quaternion `Q`.

    [procedure] (rotate-quaternion Q MATRIX)

Rotate `MATRIX` by the quaternion `Q`.

    [procedure] (ypr-rotation YAW PITCH ROLL)

Return the rotation matrix given by rotating by `ROLL` radians around the z-axis followed by `PITCH` radians around the x-axis followed by `YAW` radians around the y-axis.

    [procedure] (rotate-ypr YAW PITCH ROLL MATRIX)

Rotate `MATRIX` by `ROLL` radians around the z-axis followed by `PITCH` radians around the x-axis followed by `YAW` radians around the y-axis.

    [procedure] (2d-scaling SCALE-X SCALE-Y)

Return the matrix created by scaling the x and y axes by `SCALE-X` and `SCALE-Y`.

    [procedure] (scale-2d SCALE-X SCALE-Y MATRIX)

Scale the x and y axis of `MATRIX` by `SCALE-X` and `SCALE-Y`.

    [procedure] (3d-scaling SCALE-X SCALE-Y SCALE-Z)

Return the matrix created by scaling the x, y and z axes by `SCALE-X`, `SCALE-Y`, and `SCALE-Z`.

    [procedure] (scale-3d SCALE-X SCALE-Y SCALE-Z MATRIX)

Scale the x, y, and z axis of `MATRIX` by `SCALE-X`, `SCALE-Y`, and `SCALE-Z`.

    [procedure] (scaling SCALE)

Return the matrix created by scaling the x, y and z axes by `SCALE`.

    [procedure] (scale SCALE MATRIX)

Scale the x, y, and z axis of `MATRIX` by `SCALE`.

    [procedure] (flip-x MATRIX)

Flip (mirror) `MATRIX` along the x-axis.

    [procedure] (flip-y MATRIX)

Flip (mirror) `MATRIX` along the y-axis.

    [procedure] (flip-z MATRIX)

Flip (mirror) `MATRIX` along the z-axis.

    [procedure] (translate-rotate-scale-2d VECTOR ANGLE SCALE)

Efficiently create a matrix translated by `VECTOR`, rotated around the z-axis by `ANGLE` radians, then scaled by `SCALE`.

    [procedure] (transpose MATRIX)

Transpose `MATRIX`.

    [procedure] (inverse MATRIX)

Invert `MATRIX`.

    [procedure] (fast-inverse-transpose MATRIX)

Efficiently inverse the transpose the unscaled `MATRIX`. If `MATRIX` has been scaled, this will produce incorrect results: `inverse` then `transpose` should be used instead.

### Projection matrices
    [procedure] (ortho WIDTH HEIGHT NEAR FAR)

Create an orthographic projection matrix.

    [procedure] (ortho-viewport LEFT RIGHT NEAR FAR VIEWPORT-LEFT VIEWPORT-RIGHT VIEWPORT-BOTTOM VIEWPORT-TOP)

Create an orthographic projection matrix mapping the `LEFT`, `RIGHT`, `TOP`, `BOTTOM`, `NEAR`, `FAR` cube to a viewport of `VIEWPORT-LEFT`, `VIEWPORT-RIGHT`, `VIEWPORT-TOP`, `VIEWPORT-BOTTOM`.

    [procedure] (perspective WIDTH HEIGHT NEAR FAR ANGLE)

Create an perspective projection matrix with a field of view of `ANGLE` degrees.

    [procedure] (frustum LEFT RIGHT BOTTOM TOP NEAR FAR)

Create a perspective projection matrix defined by a frustum with a near side of `LEFT`, `RIGHT`, `TOP`, `BOTTOM`, `NEAR`, and a far side at `FAR`. If the matrix `RESULT` is given, it will be modified to contain the result. If `RESULT` is `#t`, the returned value will be an f32vector located in non-garbage-collected memory (the memory will still be freed when there are no more references to the matrix). If `RESULT` is not provided, the returned value will be an f32vector located in normal garbage collected memory.

    [procedure] (frustum-viewport LEFT RIGHT BOTTOM TOP NEAR FAR VIEWPORT-LEFT VIEWPORT-RIGHT VIEWPORT-BOTTOM VIEWPORT-TOP)

Create a perspective projection matrix mapping the `LEFT`, `RIGHT`, `TOP`, `BOTTOM`, `NEAR`, `FAR` frustum to a viewport of `VIEWPORT-LEFT`, `VIEWPORT-RIGHT`, `VIEWPORT-BOTTOM`, `VIEWPORT-TOP`.


### Camera functions
    [procedure] (look-at EYE OBJ UP)

Create a “look-at” style camera matrix. The camera is positioned at point `EYE`, pointing towards the point `OBJ`. `UP` defines the camera’s up vector.

    [procedure] (camera-inverse CAMERA)

Invert `CAMERA` in an efficient fashion. This allows the camera to be constructed in an intuitive fashion by translating and rotating before inverting in order to position the scene properly. This function is far faster than the general `inverse` function, but the matrix `CAMERA` must only be a matrix representing a translation and a rotation (no scaling).

### Vector operations
    [procedure] (make-point X Y Z)
    [procedure] (point-x POINT)
    [procedure] (point-y POINT)
    [procedure] (point-z POINT)

Vector constructor, getters, and setters.

    [procedure] (v+ A B)
 
Return the result of the addition of vectors `A` and `B`.

    [procedure] (v- A B)
 
Return the result of the subtraction of vector `B` from `A`.

    [procedure] (v* V S)
 
Return the result of the multiplication of vector `A` with scalar `S`.

    [procedure] (cross-product A B)
 
Return the result of the cross product between the vectors `A` and `B`.

    [procedure] (dot-product A B)

Return the result of the dot product between the vectors `A` and `B`.

    [procedure] (vector-magnitude V)

Return the magnitude of vector `V`.

    [procedure] (normalize V)

Normalize the vector `V`.

    [procedure] (m*vector MATRIX VECTOR)

Multiply `VECTOR` by `MATRIX`.

    [procedure] (lerp A B T)

Linear interpolation between the points `A` and `B` with the interpolation parameter `T` which must be between 0 and 1. If the vector `RESULT` is given, it will be modified to contain the result.


### Quaternion operations
Quaternions are expected to be normalized before they are used in certain functions (`quaternion-normalize` may be used to do so). All the provided functions that create quaternions, create unit quaternions. 

The order of quaternion cross-multiplication is the inverse of the “standard” order, so a quaternion that has undergone a series or rotations will represent the same rotation as a marix that has gone through the same series, in the same order.

    [procedure] (make-quaternion X Y Z W)
    [procedure] (quaternion-x POINT)
    [procedure] (quaternion-y POINT)
    [procedure] (quaternion-z POINT)
    [procedure] (quaternion-w POINT)

Quaternion constructor, getters, and setters.

    [procedure] (quaternion-normalize Q)

Destructively normalize the quaternion `Q`.

    [procedure] (quaternion-inverse Q)

Return the inverse of the unit quaternion `Q`.

    [procedure] (quaternion-cross-product A B)

Return the cross-product of the quaternions `A` and `B`.

    [procedure] (quaternion-dot-product A B)

Return the dot-product of the quaternions `A` and `B`.

    [procedure] (q+ A B)
 
Return the result of the addition of quaternions `A` and `B`.

    [procedure] (q- A B)

Return the result of the subtraction of quaternions `B` from `A`.

    [procedure] (q* V S)

Return the result of the multiplication of quaternion `A` with scalar `S`.

    [procedure] (quaternion-axis-angle-rotation AXIS ANGLE)

Return the quaternion corresponding to a rotation of `ANGLE` radians around the vector `AXIS`.

    [procedure] (quaternion-rotate-axis-angle AXIS ANGLE Q)

Rotate the quaternion `Q` by a rotation of `ANGLE` radians around the vector `AXIS`.

    [procedure] (quaternion-x-rotation ANGLE)

Return the quaternion corresponding to a rotation of `ANGLE` radians around the x-axis.

    [procedure] (quaternion-rotate-x ANGLE Q)

Rotate the quaternion `Q` by a rotation of `ANGLE` radians around the x-axis.

    [procedure] (quaternion-y-rotation ANGLE)

Return the quaternion corresponding to a rotation of `ANGLE` radians around the y-axis.

    [procedure] (quaternion-rotate-y ANGLE Q)

Rotate the quaternion `Q` by a rotation of `ANGLE` radians around the y-axis.

    [procedure] (quaternion-z-rotation ANGLE)

Return the quaternion corresponding to a rotation of `ANGLE` radians around the z-axis.

    [procedure] (quaternion-rotate-z ANGLE Q)

Rotate the quaternion `Q` by a rotation of `ANGLE` radians around the z-axis.

    [procedure] (quaternion-ypr-rotation YAW PITCH ROLL)

Return the quaternion corresponding to a rotation of `ROLL` radians around the z-axis followed by `PITCH` radians around the x-axis followed by `YAW` radians around the y-axis.

    [procedure] (quaternion-rotate-ypr YAW PITCH ROLL Q)

Rotate the quaternion `Q` by `ROLL` radians around the z-axis followed by `PITCH` radians around the x-axis followed by `YAW` radians around the y-axis.

    [procedure] (quaternion-rotate-point Q P)

Destructively rotate the point `P` by the unit quaternion `Q`.

    [procedure] (slerp A B T)

Spherical linear interpolation between the quaternions `A` and `B` with the interpolation parameter `T` which must be between 0 and 1.


### Angle operations
    [procedure] (degrees->radians ANGLE)

Converts `ANGLE` from degrees to radians.

    [procedure] (radians->degrees ANGLE)

Converts `ANGLE` from radians to degrees.

## Example

``` Scheme
(import (gl-math))

(define projection-matrix
  (perspective 640 480 0.1 100 70))

(define view-matrix
  (look-at (make-point 1 0 3)
           (make-point 0 0 0)
           (make-point 0 1 0)))

(define model-matrix (mat4-identity))

(pretty-print (m* projection-matrix
                (m* view-matrix model-matrix)))
```

## Version history
### Version 0.8.0
27 June 2020

* Forked from main.

### Version 0.8.0
7 August 2014

* Add `m*`, `m+`, `m-`

### Version 0.7.0
* Add `frustum`, `frustum-viewport`

### Version 0.6.0
8 October 2014

* Add `fast-inverse-transpose`

### Version 0.5.2
10 September 2014

* `m*vector-array!`: Stride is given in bytes when vector is a pointer

**Version 0.5.0**

2 September 2014

* Many new vector and quaternion functions
* Functions that previously accepted vectors as individual floats, now accept them as f32vectors

### Version 0.4.1
30 August 2014

* Fix `m*vector-array!`

**Version 0.4.0**

27 July 2014

* Add `copy-mat4`

### Version 0.3.2
21 July 2014

* Allow pointer to array of vectors to be passed to `m*vector-array!`
* Fix error forms

**Version 0.3.1**

23 June 2014

* Matrix vector multiplication

### Version 0.2.0
21 June 2014

* Each transformation function now has two variants: one that initializes a matrix, and one that operates on a matrix
* Provide `pi/2`
* Provide quaternion and YPR rotation
* Remove unhelpful composite operations
* Fix optional arguments for matrix operations
* Fix a bug in `look-at`

### Version 0.1.0
17 June 2014

* Initial release

## Source repository
Source available on [GitHub](https://github.com/swatson555/gl-math).
Source available on [GitHub](https://github.com/AlexCharlton/gl-math).

Bug reports and patches welcome! Bugs can be reported via GitHub.

## Author
Alex Charlton, Steven Watson

## Licence
BSD
