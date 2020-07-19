(library (gl-math)
  (export degrees->radians
	  radians->degrees

	  make-point
	  point-x
	  point-y
	  point-z
	  m*vector
	  cross-product
	  v+
	  v-
	  v*
	  point-magnitude
	  normalize
	  dot-product
	  lerp

	  make-quaternion
	  quaternion-x
	  quaternion-y
	  quaternion-z
	  quaternion-w
	  quaternion-magnitude
	  quaternion-normalize
	  quaternion-inverse
	  quaternion-dot-product
	  quaternion-cross-product
	  q+
	  q-
	  q*
	  m*quaternion
	  quaternion-rotate-point
	  quaternion-axis-angle-rotation
	  quaternion-rotate-axis-angle
	  quaternion-x-rotation
	  quaternion-rotate-x
	  quaternion-y-rotation
	  quaternion-rotate-y
	  quaternion-z-rotation
	  quaternion-rotate-z
	  quaternion-ypr-rotation
	  quaternion-rotate-ypr
	  slerp

	  mat4-identity
	  ortho-viewport
	  ortho
	  frustum-viewport
	  frustum
	  perspective
	  look-at
	  camera-inverse
	  m*
	  m*s
	  m+
	  m-
	  translation
	  translate
	  x-rotation
	  rotate-x
	  y-rotation
	  rotate-y
	  z-rotation
	  rotate-z
	  axis-angle-rotation
	  rotate-axis-angle
	  quaternion-rotation
	  rotate-quaternion
	  ypr-rotation
	  rotate-ypr
	  scaling-2d
	  scale-2d
	  scaling-3d
	  scale-3d
	  scaling
	  scale
	  flip-x
	  flip-y
	  flip-z
	  translate-rotate-scale-2d
	  transpose
	  inverse
	  fast-inverse-transpose)

  (import (rnrs))

  (define (degrees->radians degrees)
    (* degrees 0.0174532925))

  (define (radians->degrees radians)
    (* radians 57.2957795))


  (define (make-point x y z)
    `#(,x ,y ,z))

  (define (point-x p)
    (vector-ref p 0))

  (define (point-y p)
    (vector-ref p 1))

  (define (point-z p)
    (vector-ref p 2))

  (define (m*vector matrix vector)
    (define x (point-x vector))
    (define y (point-y vector))
    (define z (point-z vector))
    (define w (+ (* (_41 matrix) x)
		 (* (_42 matrix) y)
		 (* (_43 matrix) z)
		 (_44 matrix)))
    (make-point (* w (+ (* (_11 matrix) x)
		      (* (_12 matrix) y)
		      (* (_13 matrix) z)
		      (_14 matrix)))
		(* w (+ (* (_21 matrix) x)
		      (* (_22 matrix) y)
		      (* (_23 matrix) z)
		      (_24 matrix)))
		(* w (+ (* (_31 matrix) x)
		      (* (_32 matrix) y)
		      (* (_33 matrix) z)
		      (_34 matrix)))))

  (define (cross-product a b)
    (make-point (- (* (point-y a) (point-z b))
		   (* (point-z a) (point-y b)))
		(- (* (point-z a) (point-x b))
		   (* (point-x a) (point-z b)))
		(- (* (point-x a) (point-y b))
		   (* (point-y a) (point-x b)))))

  (define (v+ a b)
    (make-point (+ (point-x a) (point-x b))
		(+ (point-y a) (point-y b))
		(+ (point-z a) (point-z b))))

  (define (v- a b)
    (make-point (- (point-x a) (point-x b))
		(- (point-y a) (point-y b))
		(- (point-z a) (point-z b))))

  (define (v* v s)
    (make-point (* (point-x v) s)
		(* (point-y v) s)
		(* (point-z v) s)))

  (define (point-magnitude v)
    (sqrt (+ (expt (point-x v) 2)
	     (expt (point-y v) 2)
	     (expt (point-z v) 2))))

  (define (normalize v)
    (define len (point-magnitude v))
    (make-point (/ (point-x v) len)
		(/ (point-y v) len)
		(/ (point-z v) len)))

  (define (dot-product a b)
    (+ (* (point-x a)
	  (point-x b))
       (* (point-y a)
	  (point-y b))
       (* (point-z a)
	  (point-z b))))

  (define (lerp a b t)
    (define t0 (- 1 t))
    (make-point (+ (* (point-x a) t0)
		   (* (point-x b) t))
		(+ (* (point-y a) t0)
		   (* (point-y b) t))
		(+ (* (point-z a) t0)
		   (* (point-z b) t))))


  (define (make-quaternion x y z w)
    `#(,x ,y ,z ,w))

  (define (quaternion-x q)
    (vector-ref q 0))

  (define (quaternion-y q)
    (vector-ref q 1))

  (define (quaternion-z q)
    (vector-ref q 2))

  (define (quaternion-w q)
    (vector-ref q 3))

  (define (quaternion-magnitude q)
    (sqrt (+ (expt (quaternion-x q) 2)
	     (expt (quaternion-y q) 2)
	     (expt (quaternion-z q) 2)
	     (expt (quaternion-w q) 2))))

  (define (quaternion-normalize q)
    (define mag (quaternion-magnitude q))
    (make-quaternion (/ (quaternion-x q) mag)
		     (/ (quaternion-y q) mag)
		     (/ (quaternion-z q) mag)
		     (/ (quaternion-w q) mag)))

  (define (quaternion-inverse q)
    (make-quaternion (- (quaternion-x q))
		     (- (quaternion-y q))
		     (- (quaternion-z q))
		     (quaternion-w q)))

  (define (quaternion-dot-product a b)
    (+ (* (quaternion-x a) (quaternion-x b))
       (* (quaternion-y a) (quaternion-y b))
       (* (quaternion-z a) (quaternion-z b))
       (* (quaternion-w a) (quaternion-w b))))

  (define (quaternion-cross-product a b)
    (make-quaternion (- (+ (* (quaternion-w a) (quaternion-x b))
			   (* (quaternion-x a) (quaternion-w b))
			   (* (quaternion-y a) (quaternion-z b)))
			(* (quaternion-z a) (quaternion-y b)))
		     (- (+ (* (quaternion-w a) (quaternion-y b))
			   (* (quaternion-y a) (quaternion-w b))
			   (* (quaternion-z a) (quaternion-x b)))
			(* (quaternion-x a) (quaternion-z b)))
		     (- (+ (* (quaternion-w a) (quaternion-z b))
			   (* (quaternion-z a) (quaternion-w b))
			   (* (quaternion-x a) (quaternion-y b)))
			(* (quaternion-y a) (quaternion-x b)))
		     (- (* (quaternion-w a) (quaternion-w b))
			(* (quaternion-x a) (quaternion-x b))
			(* (quaternion-y a) (quaternion-y b))
			(* (quaternion-z a) (quaternion-z b)))))

  (define (q+ a b)
    (make-quaternion (+ (quaternion-x a) (quaternion-x b))
		     (+ (quaternion-y a) (quaternion-y b))
		     (+ (quaternion-z a) (quaternion-z b))
		     (+ (quaternion-w a) (quaternion-w b))))

  (define (q- a b)
    (make-quaternion (- (quaternion-x a) (quaternion-x b))
		     (- (quaternion-y a) (quaternion-y b))
		     (- (quaternion-z a) (quaternion-z b))
		     (- (quaternion-w a) (quaternion-w b))))
  (define (q* q s)
    (make-quaternion (* s (quaternion-x q))
		     (* s (quaternion-y q))
		     (* s (quaternion-z q))
		     (* s (quaternion-w q))))

  (define (quaternion-rotate-point q p)
    (define rotated
      (quaternion-cross-product
       (quaternion-cross-product
	(quaternion-inverse q)
	(make-quaternion (point-x p)
			 (point-y p)
			 (point-x p)
			 0))
       q))
    (make-point (quaternion-x rotated)
		(quaternion-y rotated)
		(quaternion-z rotated)))

  (define (quaternion-axis-angle-rotation axis angle)
    (define s (sin (/ angle 2)))
    (define c (cos (/ angle 2)))
    (make-quaternion (* s (point-x axis))
		     (* s (point-y axis))
		     (* s (point-z axis))
		     c))

  (define (quaternion-rotate-axis-angle axis angle quat)
    (quaternion-cross-product
     (quaternion-axis-angle-rotation axis angle) quat))

  (define (quaternion-x-rotation angle)
    (define axis (make-point 1 0 0))
    (quaternion-axis-angle-rotation axis angle))

  (define (quaternion-rotate-x angle quat)
    (quaternion-cross-product
     (quaternion-x-rotation angle) quat))

  (define (quaternion-y-rotation angle)
    (define axis (make-point 0 1 0))
    (quaternion-axis-angle-rotation axis angle))

  (define (quaternion-rotate-y angle quat)
    (quaternion-cross-product
     (quaternion-y-rotation angle) quat))

  (define (quaternion-z-rotation angle)
    (define axis (make-point 0 0 1))
    (quaternion-axis-angle-rotation axis angle))

  (define (quaternion-rotate-z angle quat)
    (quaternion-cross-product
     (quaternion-z-rotation angle) quat))

  (define (quaternion-ypr-rotation yaw pitch roll)
    (quaternion-rotate-z roll
     (quaternion-rotate-x pitch
      (quaternion-y-rotation yaw))))

  (define (quaternion-rotate-ypr yaw pitch roll quat)
    (quaternion-cross-product
     (quaternion-ypr-rotation yaw pitch roll) quat))

  (define (slerp a b t)
    (define cos-omega (+ (* (quaternion-x a) (quaternion-x b))
			 (* (quaternion-y a) (quaternion-y b))
			 (* (quaternion-z a) (quaternion-z b))
			 (* (quaternion-w a) (quaternion-w b))))
    (define *cos-omega
      (if (< cos-omega 0) (- cos-omega) cos-omega))
    (define *b
      (if (< cos-omega 0)
	  (make-quaternion (- (quaternion-x b))
			   (- (quaternion-y b))
			   (- (quaternion-z b))
			   (- (quaternion-w b)))
	  b))
    (if (> *cos-omega 0.9999)
	(let ((ka (- 1 t)) (kb t))
	  (make-quaternion (+ (* (quaternion-x a) ka) (* (quaternion-x *b) kb))
			   (+ (* (quaternion-y a) ka) (* (quaternion-y *b) kb))
			   (+ (* (quaternion-z a) ka) (* (quaternion-z *b) kb))
			   (+ (* (quaternion-w a) ka) (* (quaternion-w *b) kb))))
	(let* ((sin-omega (sqrt (- 1 (* *cos-omega *cos-omega))))
	       (one-over-sin-omega (/ 1 sin-omega))
	       (omega (atan sin-omega *cos-omega)))
	  (define ka (* (sin (* (- 1 t) omega)) one-over-sin-omega))
	  (define kb (* (sin (* t omega)) one-over-sin-omega))
	  (make-quaternion (+ (* (quaternion-x a) ka) (* (quaternion-x *b) kb))
			   (+ (* (quaternion-y a) ka) (* (quaternion-y *b) kb))
			   (+ (* (quaternion-z a) ka) (* (quaternion-z *b) kb))
			   (+ (* (quaternion-w a) ka) (* (quaternion-w *b) kb))))))

  (define (_11 mat)
    (vector-ref mat 0))

  (define (_21 mat)
    (vector-ref mat 1))

  (define (_31 mat)
    (vector-ref mat 2))

  (define (_41 mat)
    (vector-ref mat 3))

  (define (_12 mat)
    (vector-ref mat 4))

  (define (_22 mat)
    (vector-ref mat 5))

  (define (_32 mat)
    (vector-ref mat 6))

  (define (_42 mat)
    (vector-ref mat 7))

  (define (_13 mat)
    (vector-ref mat 8))

  (define (_23 mat)
    (vector-ref mat 9))

  (define (_33 mat)
    (vector-ref mat 10))

  (define (_43 mat)
    (vector-ref mat 11))

  (define (_14 mat)
    (vector-ref mat 12))

  (define (_24 mat)
    (vector-ref mat 13))

  (define (_34 mat)
    (vector-ref mat 14))

  (define (_44 mat)
    (vector-ref mat 15))

  (define (mat4-identity)
    '#(1 0 0 0
       0 1 0 0
       0 0 1 0
       0 0 0 1))

  (define (ortho-viewport left right bottom top near far view-left view-right view-bottom view-top)
    (define *11 (/ (- view-right view-left) (- right left)))
    (define *14 (/ (- (* right view-left) (* view-right left)) (- right left)))
    (define *22 (/ (- view-top view-bottom) (- top bottom)))
    (define *24 (/ (- (* top view-bottom) (* view-top bottom)) (- top bottom)))
    (define *33 (/ -2 (- far near)))
    (define *34 (- (/ (+ far near) (- far near))))
    (define *44 1)
    `#(,*11    0    0    0
	  0 ,*22    0    0
	  0    0 ,*33    0
       ,*14 ,*24 ,*34 ,*44))

  (define (ortho width height near far)
    (define right (/ width 2))
    (define left (- right))
    (define top (/ height 2))
    (define bottom (- top))
    (ortho-viewport left right bottom top near far -1 1 1 -1))

  (define (frustum-viewport left right bottom top near far view-left view-right view-bottom view-top)
    (define *11 (/ (* near (- view-right view-left)) (- right left)))
    (define *13 (/ (- (* view-right left) (* right view-left)) (- right left)))
    (define *22 (/ (* near (- view-top view-bottom)) (- top bottom)))
    (define *23 (/ (- (* view-top bottom) (* top view-bottom)) (- top bottom)))
    (define *33 (- (/ (+ far near) (- far near))))
    (define *34 (- (/ (* 2 far near) (- far near))))
    (define *43 -1)
    `#(,*11    0    0    0
          0 ,*22    0    0
       ,*13 ,*23 ,*33 ,*43
          0    0 ,*34    0))

  (define (frustum left right bottom top near far)
    (frustum-viewport left right bottom top near far -1 1 -1 1))

  (define (perspective width height near far fov-angle)
    (define top (* (tan (* fov-angle 0.5)) near))
    (define bottom (- top))
    (define right (* (/ width height) top))
    (define left (- right))
    (frustum left right bottom top near far))

  (define (look-at eye obj up)
    (define f (normalize (v- eye obj)))
    (define r (normalize (cross-product up f)))
    (define u (cross-product f r))
    `#(,(point-x r) ,(point-x u) ,(point-x f) 0
       ,(point-y r) ,(point-y u) ,(point-y f) 0
       ,(point-z r) ,(point-z u) ,(point-z f) 0
       ,(- (dot-product r eye)) ,(- (dot-product u eye)) ,(- (dot-product f eye)) 1))

  (define (camera-inverse camera)
    (define *14
      (- (+ (* (_11 camera) (_14 camera))
	    (* (_21 camera) (_24 camera))
	    (* (_31 camera) (_34 camera)))))
    (define *24
      (- (+ (* (_12 camera) (_14 camera))
	    (* (_22 camera) (_24 camera))
	    (* (_32 camera) (_34 camera)))))
    (define *34
      (- (+ (* (_13 camera) (_14 camera))
	    (* (_23 camera) (_24 camera))
	    (* (_33 camera) (_34 camera)))))
    `#(,(_11 camera) ,(_12 camera) ,(_13 camera) 0
       ,(_21 camera) ,(_22 camera) ,(_23 camera) 0
       ,(_31 camera) ,(_32 camera) ,(_33 camera) 0
                ,*14          ,*24          ,*34 1))


  (define (m* a b)
    `#(,(+ (* (_11 a) (_11 b)) (* (_12 a) (_21 b)) (* (_13 a) (_31 b)) (* (_14 a) (_41 b)))
       ,(+ (* (_21 a) (_11 b)) (* (_22 a) (_21 b)) (* (_23 a) (_31 b)) (* (_24 a) (_41 b)))
       ,(+ (* (_31 a) (_11 b)) (* (_32 a) (_21 b)) (* (_33 a) (_31 b)) (* (_34 a) (_41 b)))
       ,(+ (* (_41 a) (_11 b)) (* (_42 a) (_21 b)) (* (_43 a) (_31 b)) (* (_44 a) (_41 b)))
       ,(+ (* (_11 a) (_12 b)) (* (_12 a) (_22 b)) (* (_13 a) (_32 b)) (* (_14 a) (_42 b)))
       ,(+ (* (_21 a) (_12 b)) (* (_22 a) (_22 b)) (* (_23 a) (_32 b)) (* (_24 a) (_42 b)))
       ,(+ (* (_31 a) (_12 b)) (* (_32 a) (_22 b)) (* (_33 a) (_32 b)) (* (_34 a) (_42 b)))
       ,(+ (* (_41 a) (_12 b)) (* (_42 a) (_22 b)) (* (_43 a) (_32 b)) (* (_44 a) (_42 b)))
       ,(+ (* (_11 a) (_13 b)) (* (_12 a) (_23 b)) (* (_13 a) (_33 b)) (* (_14 a) (_43 b)))
       ,(+ (* (_21 a) (_13 b)) (* (_22 a) (_23 b)) (* (_23 a) (_33 b)) (* (_24 a) (_43 b)))
       ,(+ (* (_31 a) (_13 b)) (* (_32 a) (_23 b)) (* (_33 a) (_33 b)) (* (_34 a) (_43 b)))
       ,(+ (* (_41 a) (_13 b)) (* (_42 a) (_23 b)) (* (_43 a) (_33 b)) (* (_44 a) (_43 b)))
       ,(+ (* (_11 a) (_14 b)) (* (_12 a) (_24 b)) (* (_13 a) (_34 b)) (* (_14 a) (_44 b)))
       ,(+ (* (_21 a) (_14 b)) (* (_22 a) (_24 b)) (* (_23 a) (_34 b)) (* (_24 a) (_44 b)))
       ,(+ (* (_31 a) (_14 b)) (* (_32 a) (_24 b)) (* (_33 a) (_34 b)) (* (_34 a) (_44 b)))
       ,(+ (* (_41 a) (_14 b)) (* (_42 a) (_24 b)) (* (_43 a) (_34 b)) (* (_44 a) (_44 b)))))

  (define (m*s mat s)
    `#(,(* s (_11 mat)) ,(* s (_21 mat)) ,(* s (_31 mat)) ,(* s (_41 mat))
       ,(* s (_12 mat)) ,(* s (_22 mat)) ,(* s (_32 mat)) ,(* s (_42 mat))
       ,(* s (_13 mat)) ,(* s (_23 mat)) ,(* s (_33 mat)) ,(* s (_43 mat))
       ,(* s (_14 mat)) ,(* s (_24 mat)) ,(* s (_34 mat)) ,(* s (_44 mat))))

  (define (m+ a b)
    `#(,(+ (_11 a) (_11 b)) ,(+ (_21 a) (_21 b)) ,(+ (_31 a) (_31 b)) ,(+ (_41 a) (_41 b))
       ,(+ (_12 a) (_12 b)) ,(+ (_22 a) (_22 b)) ,(+ (_32 a) (_32 b)) ,(+ (_42 a) (_42 b))
       ,(+ (_13 a) (_13 a)) ,(+ (_23 a) (_23 a)) ,(+ (_33 a) (_33 a)) ,(+ (_43 a) (_43 a))
       ,(+ (_14 a) (_14 b)) ,(+ (_24 a) (_24 b)) ,(+ (_34 a) (_34 b)) ,(+ (_44 a) (_44 b))))

  (define (m- a b)
    `#(,(- (_11 a) (_11 b)) ,(- (_21 a) (_21 b)) ,(- (_31 a) (_31 b)) ,(- (_41 a) (_41 b))
       ,(- (_12 a) (_12 b)) ,(- (_22 a) (_22 b)) ,(- (_32 a) (_32 b)) ,(- (_42 a) (_42 b))
       ,(- (_13 a) (_13 a)) ,(- (_23 a) (_23 a)) ,(- (_33 a) (_33 a)) ,(- (_43 a) (_43 a))
       ,(- (_14 a) (_14 b)) ,(- (_24 a) (_24 b)) ,(- (_34 a) (_34 b)) ,(- (_44 a) (_44 b))))

  (define (translation p)
    (define x (point-x p))
    (define y (point-y p))
    (define z (point-z p))
    `#( 1  0  0 0
        0  1  0 0
        0  0  1 0
       ,x ,y ,z 1))

  (define (translate p matrix)
    (m* (translation p) matrix))

  (define (x-rotation rotation)
    (define c (cos rotation))
    (define s (sin rotation))
    `#(1      0  0 0
       0     ,c ,s 0
       0 ,(- s) ,c 0
       0      0  0 1))

  (define (rotate-x rotation mat)
    (m* (x-rotation rotation) mat))

  (define (y-rotation rotation)
    (define c (cos rotation))
    (define s (sin rotation))
    `#(,c 0 ,(- s) 0
        0 1      0 0
       ,s 0     ,c 0
        0 0      0 1))

  (define (rotate-y rotation mat)
    (m* (y-rotation rotation) mat))

  (define (z-rotation rotation)
    (define c (cos rotation))
    (define s (sin rotation))
    `#(    ,c ,s 0 0
       ,(- s) ,c 0 0
            0  0 1 0
            0  0 0 1))

  (define (rotate-z rotation mat)
    (m* (z-rotation rotation) mat))

  (define (axis-angle-rotation axis angle)
    (define c (cos angle))
    (define s (sin angle))
    (define C (- 1 c))
    (define a (normalize axis))
    (define xx (* (point-x a) (point-x a)))
    (define xy (* (point-x a) (point-y a)))
    (define xz (* (point-x a) (point-z a)))
    (define yy (* (point-y a) (point-y a)))
    (define yz (* (point-y a) (point-z a)))
    (define zz (* (point-z a) (point-z a)))
    (define xs (* (point-x a) s))
    (define ys (* (point-y a) s))
    (define zs (* (point-z a) s))
    `#(,(+ (* xx C) c)  ,(+ (* xy C) zs) ,(- (* xz C) ys) 0
       ,(- (* xy C) zs) ,(+ (* yy C) c)  ,(+ (* yz C) xs) 0
       ,(+ (* xz C) ys) ,(- (* yz C) xs) ,(+ (* zz C) c)  0
       0 0 0 1))

  (define (rotate-axis-angle axis angle matrix)
    (m* (axis-angle-rotation axis angle) matrix))

  (define (quaternion-rotation q)
    (define xx (* (quaternion-x q) (quaternion-x q)))
    (define xy (* (quaternion-x q) (quaternion-y q)))
    (define xz (* (quaternion-x q) (quaternion-z q)))
    (define xw (* (quaternion-x q) (quaternion-w q)))
    (define yy (* (quaternion-y q) (quaternion-y q)))
    (define yz (* (quaternion-y q) (quaternion-z q)))
    (define yw (* (quaternion-y q) (quaternion-w q)))
    (define zz (* (quaternion-z q) (quaternion-z q)))
    (define zw (* (quaternion-z q) (quaternion-w q)))
    `#(,(- 1 (* 2 (+ zz yy))) ,(* 2 (+ xy zw))       ,(* 2 (- xz yw))       0
       ,(* 2 (- xy zw))       ,(- 1 (* 2 (+ xx zz))) ,(* 2 (+ yz xw))       0
       ,(* 2 (+ xz yw))       ,(* 2 (* (- yz xw)))   ,(- 1 (* 2 (+ xx yy))) 0
       0 0 0 1))

  (define (rotate-quaternion q matrix)
    (m* (quaternion-rotation q) matrix))

  (define (ypr-rotation yaw pitch roll)
    (define sy (sin yaw))
    (define cy (cos yaw))
    (define sp (sin pitch))
    (define cp (cos pitch))
    (define sr (sin roll))
    (define cr (cos roll))
    `#(,(+ (* sp sy sr) (* cy cr)) ,(* cp sr) ,(- (* sp cy sr) (* sy cr)) 0
       ,(- (* sp sy cr) (* cy sr)) ,(* cp cr) ,(+ (* sp cy cr) (* sy sr)) 0
       ,(* cp sy)                  ,(- sp)    ,(* cp cy)                  0
       0 0 0 1))

  (define (rotate-ypr yaw pitch roll matrix)
    (m* (ypr-rotation yaw pitch roll) matrix))

  (define (scaling-2d scale-x scale-y)
    `#(,scale-x        0 0 0
              0 ,scale-y 0 0
	      0        0 1 0
	      0        0 0 1))

  (define (scale-2d scale-x scale-y matrix)
    (m* (scaling-2d scale-x scale-y) matrix))

  (define (scaling-3d scale-x scale-y scale-z)
    `#(,scale-x        0        0 0
              0 ,scale-y        0 0
	      0        0 ,scale-z 0
	      0        0        0 1))

  (define (scale-3d scale-x scale-y scale-z matrix)
    (m* (scaling-3d scale-x scale-y scale-z) matrix))

  (define (scaling factor)
    `#(,factor       0       0 0
             0 ,factor       0 0
	     0       0 ,factor 0
	     0       0       0 1))

  (define (scale factor matrix)
    (m* (scaling factor) matrix))

  (define (flip-x matrix)
    `#(,(_11 matrix) ,(_21 matrix)     ,(_31 matrix) ,(_41 matrix)
       ,(_12 matrix) ,(- (_22 matrix)) ,(_32 matrix) ,(_42 matrix)
       ,(_13 matrix) ,(_23 matrix)     ,(_33 matrix) ,(_43 matrix)
       ,(_14 matrix) ,(_24 matrix)     ,(_34 matrix) ,(_44 matrix)))

  (define (flip-y matrix)
    `#(,(- (_11 matrix)) ,(_21 matrix) ,(_31 matrix) ,(_41 matrix)
       ,(_12 matrix)     ,(_22 matrix) ,(_32 matrix) ,(_42 matrix)
       ,(_13 matrix)     ,(_23 matrix) ,(_33 matrix) ,(_43 matrix)
       ,(_14 matrix)     ,(_24 matrix) ,(_34 matrix) ,(_44 matrix)))

  (define (flip-z matrix)
    `#(,(_11 matrix) ,(_21 matrix) ,(_31 matrix)     ,(_41 matrix)
       ,(_12 matrix) ,(_22 matrix) ,(_32 matrix)     ,(_42 matrix)
       ,(_13 matrix) ,(_23 matrix) ,(- (_33 matrix)) ,(_43 matrix)
       ,(_14 matrix) ,(_24 matrix) ,(_34 matrix)     ,(_44 matrix)))

  (define (translate-rotate-scale-2d p angle scale)
    (define c (* scale (cos angle)))
    (define s (* scale (sin angle)))
    `#(,c     ,s  0 0
       ,(- s) ,c  0 0
       0      0   1 0
       ,(point-x p) ,(point-y p) ,(point-z p) 1))

  (define (transpose matrix)
    `#(,(_11 matrix) ,(_12 matrix) ,(_13 matrix) ,(_14 matrix)
       ,(_21 matrix) ,(_22 matrix) ,(_23 matrix) ,(_24 matrix)
       ,(_31 matrix) ,(_32 matrix) ,(_33 matrix) ,(_34 matrix)
       ,(_41 matrix) ,(_42 matrix) ,(_43 matrix) ,(_44 matrix)))

  (define (inverse matrix)
    (define *11 (- (+ (- (* (_22 matrix) (_33 matrix) (_44 matrix))
			 (* (_22 matrix) (_34 matrix) (_43 matrix))
			 (* (_32 matrix) (_23 matrix) (_44 matrix)))
		      (* (_32 matrix) (_24 matrix) (_43 matrix))
		      (* (_42 matrix) (_23 matrix) (_34 matrix)))
		   (* (_42 matrix) (_24 matrix) (_33 matrix))))
    (define *21 (+ (- (+ (- (* (_21 matrix) (_33 matrix) (_44 matrix)))
			 (* (_21 matrix) (_34 matrix) (_43 matrix))
			 (* (_31 matrix) (_23 matrix) (_44 matrix)))
		      (* (_31 matrix) (_24 matrix) (_43 matrix))
		      (* (_41 matrix) (_23 matrix) (_34 matrix)))
		   (* (_41 matrix) (_24 matrix) (_33 matrix))))
    (define *31 (- (+ (- (* (_21 matrix) (_32 matrix) (_44 matrix))
			 (* (_21 matrix) (_34 matrix) (_42 matrix))
			 (* (_31 matrix) (_22 matrix) (_44 matrix)))
		      (* (_31 matrix) (_24 matrix) (_42 matrix))
		      (* (_41 matrix) (_22 matrix) (_34 matrix)))
		   (* (_41 matrix) (_24 matrix) (_32 matrix))))
    (define *41 (+ (- (+ (- (* (_21 matrix) (_32 matrix) (_43 matrix)))
			 (* (_21 matrix) (_33 matrix) (_42 matrix))
			 (* (_31 matrix) (_22 matrix) (_43 matrix)))
		      (* (_31 matrix) (_23 matrix) (_42 matrix))
		      (* (_41 matrix) (_22 matrix) (_33 matrix)))
		   (* (_41 matrix) (_23 matrix) (_32 matrix))))
    (define *12 (+ (- (+ (- (* (_12 matrix) (_33 matrix) (_44 matrix)))
			 (* (_12 matrix) (_34 matrix) (_43 matrix))
			 (* (_32 matrix) (_13 matrix) (_44 matrix)))
		      (* (_32 matrix) (_14 matrix) (_43 matrix))
		      (* (_42 matrix) (_13 matrix) (_34 matrix)))
		   (* (_42 matrix) (_14 matrix) (_33 matrix))))
    (define *22 (- (+ (- (* (_11 matrix) (_33 matrix) (_44 matrix))
			 (* (_11 matrix) (_34 matrix) (_43 matrix))
			 (* (_31 matrix) (_13 matrix) (_44 matrix)))
		      (* (_31 matrix) (_14 matrix) (_43 matrix))
		      (* (_41 matrix) (_13 matrix) (_34 matrix)))
		   (* (_41 matrix) (_14 matrix) (_33 matrix))))
    (define *32 (+ (- (+ (- (* (_11 matrix) (_32 matrix) (_44 matrix)))
			 (* (_11 matrix) (_34 matrix) (_42 matrix))
			 (* (_31 matrix) (_12 matrix) (_44 matrix)))
		      (* (_31 matrix) (_14 matrix) (_42 matrix))
		      (* (_41 matrix) (_12 matrix) (_34 matrix)))
		   (* (_41 matrix) (_14 matrix) (_32 matrix))))
    (define *42 (- (+ (- (* (_11 matrix) (_32 matrix) (_43 matrix))
			 (* (_11 matrix) (_33 matrix) (_42 matrix))
			 (* (_31 matrix) (_12 matrix) (_43 matrix)))
		      (* (_31 matrix) (_13 matrix) (_42 matrix))
		      (* (_41 matrix) (_12 matrix) (_33 matrix)))
		   (* (_41 matrix) (_13 matrix) (_32 matrix))))
    (define *13 (- (+ (- (* (_12 matrix) (_23 matrix) (_44 matrix))
			 (* (_12 matrix) (_24 matrix) (_43 matrix))
			 (* (_22 matrix) (_13 matrix) (_44 matrix)))
		      (* (_22 matrix) (_14 matrix) (_43 matrix))
		      (* (_42 matrix) (_13 matrix) (_24 matrix)))
		   (* (_42 matrix) (_14 matrix) (_23 matrix))))
    (define *23 (+ (- (+ (- (* (_11 matrix) (_23 matrix) (_44 matrix)))
			 (* (_11 matrix) (_24 matrix) (_43 matrix))
			 (* (_21 matrix) (_13 matrix) (_44 matrix)))
		      (* (_21 matrix) (_14 matrix) (_43 matrix))
		      (* (_41 matrix) (_13 matrix) (_24 matrix)))
		   (* (_41 matrix) (_14 matrix) (_23 matrix))))
    (define *33 (- (+ (- (* (_11 matrix) (_22 matrix) (_44 matrix))
			 (* (_11 matrix) (_24 matrix) (_42 matrix))
			 (* (_21 matrix) (_12 matrix) (_44 matrix)))
		      (* (_21 matrix) (_14 matrix) (_42 matrix))
		      (* (_41 matrix) (_12 matrix) (_24 matrix)))
		   (* (_41 matrix) (_14 matrix) (_22 matrix))))
    (define *43 (+ (- (+ (- (* (_11 matrix) (_22 matrix) (_43 matrix)))
			 (* (_11 matrix) (_23 matrix) (_42 matrix))
			 (* (_21 matrix) (_12 matrix) (_43 matrix)))
		      (* (_21 matrix) (_13 matrix) (_42 matrix))
		      (* (_41 matrix) (_12 matrix) (_23 matrix)))
		   (* (_41 matrix) (_13 matrix) (_22 matrix))))
    (define *14 (+ (- (+ (- (* (_12 matrix) (_23 matrix) (_34 matrix)))
			 (* (_12 matrix) (_24 matrix) (_33 matrix))
			 (* (_22 matrix) (_13 matrix) (_34 matrix)))
		      (* (_22 matrix) (_14 matrix) (_33 matrix))
		      (* (_32 matrix) (_13 matrix) (_24 matrix)))
		   (* (_32 matrix) (_14 matrix) (_23 matrix))))
    (define *24 (- (+ (- (* (_11 matrix) (_23 matrix) (_34 matrix))
			 (* (_11 matrix) (_24 matrix) (_33 matrix))
			 (* (_21 matrix) (_13 matrix) (_34 matrix)))
		      (* (_21 matrix) (_14 matrix) (_33 matrix))
		      (* (_31 matrix) (_13 matrix) (_24 matrix)))
		   (* (_31 matrix) (_14 matrix) (_23 matrix))))
    (define *34 (+ (- (+ (- (* (_11 matrix) (_22 matrix) (_34 matrix)))
			 (* (_11 matrix) (_24 matrix) (_32 matrix))
			 (* (_21 matrix) (_12 matrix) (_34 matrix)))
		      (* (_21 matrix) (_14 matrix) (_32 matrix))
		      (* (_31 matrix) (_12 matrix) (_24 matrix)))
		   (* (_31 matrix) (_14 matrix) (_22 matrix))))
    (define *44 (- (+ (- (* (_11 matrix) (_22 matrix) (_33 matrix))
			 (* (_11 matrix) (_23 matrix) (_32 matrix))
			 (* (_21 matrix) (_12 matrix) (_33 matrix)))
		      (* (_21 matrix) (_13 matrix) (_32 matrix))
		      (* (_31 matrix) (_12 matrix) (_23 matrix)))
		   (* (_31 matrix) (_13 matrix) (_22 matrix))))

    (define det (/ 1 (+ (* (_11 matrix) *11) (* (_12 matrix) *21) (* (_13 matrix) *31) (* (_14 matrix) *41))))

    `#(,(* *11 det) ,(* *21 det) ,(* *31 det) ,(* *41 det)
       ,(* *12 det) ,(* *22 det) ,(* *32 det) ,(* *42 det)
       ,(* *13 det) ,(* *23 det) ,(* *33 det) ,(* *43 det)
       ,(* *14 det) ,(* *24 det) ,(* *34 det) ,(* *44 det)))

  (define (fast-inverse-transpose matrix)
    (define *41 (- (+ (* (_11 matrix) (_14 matrix)) (* (_21 matrix) (_24 matrix)) (* (_31 matrix) (_34 matrix)))))
    (define *42 (- (+ (* (_12 matrix) (_14 matrix)) (* (_22 matrix) (_24 matrix)) (* (_32 matrix) (_34 matrix)))))
    (define *43 (- (+ (* (_13 matrix) (_14 matrix)) (* (_23 matrix) (_24 matrix)) (* (_33 matrix) (_34 matrix)))))
    `#(,(_11 matrix) ,(_21 matrix) ,(_31 matrix) ,*41
       ,(_12 matrix) ,(_22 matrix) ,(_32 matrix) ,*42
       ,(_13 matrix) ,(_23 matrix) ,(_33 matrix) ,*43
       0 0 0 1)))
