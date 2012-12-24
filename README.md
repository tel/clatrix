# Clatrix matrix

A stupid name for a smart matrix library, because who doesn't love
smart matrices? Being implemented as a data type around the native
BLAS hooks of [jblas](http://github.com/mikiobraun/jblas) gives it
speed. Being implemented as a Clojure sequence makes it clever.

## Usage

For now, you can read the 
[Marginalia documentation](http://quantisan.github.com/clatrix) or take a look
at a few examples below. (Note: To be updated!)

```Clojure
(in-ns 'clatrix.core)

(let [A (matrix (repeat 5 (range 10)))] (pp A))

; A 5x10 matrix
; -------------
; 0.00e+00  1.00e+00  2.00e+00  .  7.00e+00  8.00e+00  9.00e+00 
; 0.00e+00  1.00e+00  2.00e+00  .  7.00e+00  8.00e+00  9.00e+00 
; 0.00e+00  1.00e+00  2.00e+00  .  7.00e+00  8.00e+00  9.00e+00 
; 0.00e+00  1.00e+00  2.00e+00  .  7.00e+00  8.00e+00  9.00e+00 
; 0.00e+00  1.00e+00  2.00e+00  .  7.00e+00  8.00e+00  9.00e+00 


(pp (from-indices 5 5 *))

; A 5x5 matrix
; ------------
; 0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00 
; 0.00e+00  1.00e+00  2.00e+00  3.00e+00  4.00e+00 
; 0.00e+00  2.00e+00  4.00e+00  6.00e+00  8.00e+00 
; 0.00e+00  3.00e+00  6.00e+00  9.00e+00  1.20e+01 
; 0.00e+00  4.00e+00  8.00e+00  1.20e+01  1.60e+01 

(time (solve (rand 1000 1000) (rand 1000)))

; "Elapsed time: 158.216 msecs"
; #<Matrix [1000 1]>

(time (rank (rand 1000 1000)))

; "Elapsed time: 13469.51 msecs"
; 1000

(let [A (rand 10 14)
      B (* A (t A))     ; B is symmetric
      lu (lu B)
      P (rspectral 10)  ; `respectral` makes positive definite matrices
      G (cholesky P)]   ; so we can get their square root
  (and (= A (c/t (c/t A)))
       (= A (hstack (cols A)))
       (= B (* (:p lu) (:l lu) (:u lu)))
       (= P (* (t G) G))))

; true
```

-----
Copyright © 2012 Paul Lam, Joseph Abrahamson (see LICENSE.txt)