# Clatrix matrix

A stupid name for a smart matrix library, because Clojure needs smart
matrices. Being implemented as a thin wrapper around the native BLAS
hooks of [jblas](http://github.com/mikiobraun/jblas) gives it
speed. Being implemented in Clojure makes it clever.

## Usage

For now, you can read the 
[Marginalia documentation](http://tel.github.com/clatrix) or take a look
at a few examples below.

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
      B (* A (t A))
      P (rspectral 10)
      G (cholesky P)
      lu (lu B)]
  (and (= A (c/t (c/t A)))
       (= A (hstack (cols A)))
       (= B (* (:p lu) (:l lu) (:u lu)))
       (= P (* (t G) G))))

; true
```

## License

Copyright © 2012 Joseph Abrahamson

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.