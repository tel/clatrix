(ns clatrix.core-test
  (:use expectations)
  (:require [clatrix.core :as c])
  (:import [clatrix.core Matrix]
           [java.io StringReader PushbackReader]))

(defn read* [str]
  (read (PushbackReader. (StringReader. str))))

(defn truthy? [v]
  (not (or (false? v) (nil? v))))

(def n 10)
(def m 15)
(def A (c/rnorm n m))
(def B (c/* (c/t A) A))
(def E (empty A))
(def S (c/rnorm n n))
(def p (c/rnorm m))
(def q (c/rnorm n))
(def F (c/matrix (partition m (range (* n m)))))
(def R (c/t (c/matrix (range m))))
(def C (c/column (range n)))
(def M (c/matrix [[1 2 3] [4 5 6]]))
(def ridx (range n))
(def cidx (range m))

;; properties of A/S
(expect Matrix A)
(expect Matrix B)
(expect Matrix E)
(expect Matrix S)
(expect Matrix p)
(expect Matrix q)
(expect Matrix F)
(expect Matrix R)
(expect Matrix C)
(expect Matrix M)

(given A
       (expect c/size [n m]
               c/matrix? true
               c/square? false
               c/column? false
               c/vector? false
               c/row?    false))

(given E
       (expect c/size [0 0]
               c/matrix? true
               c/square? true
               c/column? false
               c/vector? false
               c/row? false))

(given S
       (expect c/size [n n]
               c/matrix? true
               c/square? true
               c/column? false
               c/vector? false
               c/row?    false))

(given p
       (expect c/size [m 1]
               c/matrix? true
               c/square? false
               c/column? true
               c/vector? true
               c/row?    false))

(given q
       (expect c/size [n 1]
               c/matrix? true
               c/square? false
               c/column? true
               c/vector? true
               c/row?    false))

(given F
       (expect c/size [n m]
               c/matrix? true
               c/square? false
               c/column? false
               c/vector? false
               c/row? false
               count n))
(given R
       (expect c/size [1 15]
               c/matrix? true
               c/square? false
               c/column? false
               c/vector? true
               c/row? true
               count m))
(given C
       (expect c/size [10 1]
               c/matrix? true
               c/square? false
               c/column? true
               c/vector? true
               c/row? false
               count n))

(let [z (rand)]
  (c/set A 0 0 z)
  (expect z (c/get A 0 0)))

(expect 3.0 (c/get F 0 3))
(expect (c/matrix [[17 18]]) (c/get F 1 [2 3]))
(expect (c/matrix [[2] [17] [32]]) (c/get F [0 1 2] 2))
(expect (c/matrix [[3 5 7]
                   [108 110 112]
                   [138 140 142]])
        (c/get F [0 7 9] [3 5 7]))

(expect nil (meta M))

;; equality
(expect R (range m))
(expect (double m) (c/matrix [m]))

;; clojure sequence methods
(expect (c/matrix [(range 1 m)]) (rest R))
(expect (c/matrix (range 1 n)) (rest C))
(expect (c/matrix (vector (range m))) (first F))
(expect 0.0 (ffirst F))
(expect (c/matrix [[4 5 6]]) (rest M))
(expect (c/matrix []) (rest (rest M)))
(expect nil (next (next M)))
(expect (c/column (range 1 n)) (rest C))
(expect (c/matrix (partition m (range m (* n m)))) (rest F))

(expect (conj M [10 11 12]) (c/matrix [[1 2 3] [4 5 6] [10 11 12]]))
(expect (conj M [[10 11 12]]) (c/matrix [[1 2 3] [4 5 6] [10 11 12]]))
(expect (conj M M) (c/matrix [[1 2 3] [4 5 6]
                              [1 2 3] [4 5 6]]))
(expect [(c/matrix [[1 2 3]]) (c/matrix [[4 5 6]])] (map identity M))

(expect (c/matrix [(range 1 m)]) (rest R))
(expect (double (reduce + (range m))) (reduce + R))

(expect `(c/matrix ~(map #(map double %) [[1 2] [3 4]]))
        (read* (str (c/matrix [[1 2] [3 4]]))))

(expect 7.0 (nth C 7))
(expect (c/matrix [(range (* m 4) (* m (inc 4)))]) (nth F 4))
(expect [(c/matrix [4 19 34 49 64 79 94 109 124 139]) 
         (c/matrix [7 22 37 52 67 82 97 112 127 142])
         (c/matrix [6 21 36 51 66 81 96 111 126 141])]
        (c/cols F [4 7 6]))

;; properties of id
(given (c/id n)
       (expect c/size [n n]
               c/trace (double n)))

;; conversion from Clojure types is invertible
(expect A (c/matrix (c/dense A)))

;; `as-vec` knows about columns and rows
(expect (c/as-vec p) (flatten (c/as-vec p)))
(expect (c/as-vec q) (flatten (c/as-vec q)))
(expect false? (= (c/as-vec A) (flatten (c/as-vec A))))
(expect false? (= (c/as-vec S) (flatten (c/as-vec S))))

(expect (map double (range 10)) (c/as-vec (c/column (range 10))))

;; diagonal structure becomes 2-parity involutive
(expect (c/diag A) (c/diag (c/diag (c/diag A))))

;; other diagonal invariants
(expect (c/diag A) (c/diag (c/t A)))
(expect (c/trace S) (c/trace (c/t S)))

;; structure algebraic constraints
(expect A (c/t (c/t A)))
(expect A (apply c/hstack (c/cols A)))
(expect A (apply c/vstack (c/rows A)))

;; constants
(expect (double (* n m)) (reduce + (map (partial reduce +)
                                        (c/dense (c/constant n m 1)))))
(expect (double (* n m 5)) (reduce + (map (partial reduce +)
                                          (c/dense (c/constant n m 5)))))
(expect (double (* n m)) (reduce + (map (partial reduce +)
                                        (c/dense (c/ones n m)))))
(expect (double 0) (reduce + (map (partial reduce +)
                                  (c/dense (c/zeros n m)))))

(expect (double n) (reduce + (map (partial reduce +) (c/dense (c/id n)))))
(expect (double (* n 5)) (reduce + (map (partial reduce +)
                                        (c/dense (c/* 5 (c/id n))))))
(expect (double (* n 5)) (reduce + (map (partial reduce +)
                                        (c/dense (c/map (partial * 5) (c/id n))))))

;; norm and normalize
(expect (double m)
        (reduce + (map c/norm (c/cols (c/normalize A)))))

;; permutions
(expect A (c/permute A :r ridx))
(expect A (c/permute A :c cidx))
(expect A (c/permute (c/permute A :r (reverse ridx))
                     :r (reverse ridx)))
(expect A (c/permute (c/permute A :c (reverse cidx))
                     :c (reverse cidx)))

;; block matrices
(expect (c/id 30)
        (let [I (c/id 10)]
          (c/block [[I . .]
                    [_ I _]
                    [* * I]])))
(let [[c1 c2 c3] (c/cols A [0 1 2])
      Ac (c/hstack c1 c2 c3)
      [r1 r2 r3] (c/rows A [0 1 2])
      Ar (c/vstack r1 r2 r3)]
  (expect Ac (c/block [[c1 c2 c3]]))
  (expect Ar (c/block [[r1]
                       [r2]
                       [r3]])))
;; linear algebra
(expect A (c/- (c/+ A A) A))
(expect [m m] (c/size (c/* (c/t A) A)))
(expect [n n] (c/size (c/* A (c/t A))))
(expect A (c/* (c/id n) A))
(expect A (c/* A (c/id m)))
(expect A (c/* (c/id n) A (c/id m)))

(expect (c/constant n m 0) (c/+ A (c/- A)))

;; arithmetic
(expect (c/matrix [[2 4 6] [8 10 12]]) (c/mult M 2))
(expect (c/matrix [[1 4 9] [16 25 36]]) (c/mult M M))
(expect (c/matrix [[0.5 1.0 1.5] [2.0 2.5 3.0]]) (c/div M 2))
(expect (c/constant 2 3 1.0) (c/div M M))

;; LU decomposition
(let [lu (c/lu B)]
  (expect B (c/* (:p lu) (:l lu) (:u lu))))

;; Eigen decomposition
(let [Ss   (c/symmetric (c/* S (c/t S)))
      seig (c/eigen Ss)
      Co   (c/matrix [[0 2] [-2 0]]) ; has complex eigenvalues/vectors
      ceig (c/eigen Co)]
  (expect nil? (:ivalues  seig))
  (expect nil? (:ivectors seig))
  ;; Symmetric Eigenreconstruction is easy
  (expect Ss (c/* (:vectors seig)
                  (c/diag (:values seig))
                  (c/t (:vectors seig))))
  (expect truthy? (:ivalues  ceig))
  (expect truthy? (:ivectors ceig))

  ;; Asymmetric eigenreconstruction
  ;;
  ;; To see this formula, consider distributing matrix
  ;; multiplication over the rectangular forms of the complex
  ;; eigenvector (V) and eigenvalue (L) matrices
  ;;
  ;; = Q D Q*
  ;; = (V + iV I) (L + iL I) (Vt - iVt I)
  ;; = I (iV iL iVt + iV  L  Vt +  V iL  Vt -  V  L iVt)
  ;; +   ( V  L  Vt +  V iL iVt + iV  L iVt - iV iL  Vt)
  ;;
  ;; where `iA` is the imaginary part of some matrix with real part
  ;; `A` and I the imaginary unit.
  ;;
  ;; The complex part must vanish (since we cannot represent complex
  ;; matrices to take their eigensystems in the first place) so we
  ;; take only the second sum, which, notably, reduces to the normal
  ;; form if the imaginary parts vanish.
  (expect Co (c/+ (c/* (:vectors ceig)
                       (c/diag (:values ceig))
                       (c/t (:vectors ceig)))
                  (c/* (:vectors ceig)
                       (c/diag (:ivalues ceig))
                       (c/t (:ivectors ceig)))
                  (c/* (:ivectors ceig)
                       (c/diag (:values ceig))
                       (c/t (:ivectors ceig)))
                  (c/* -1
                       (:ivectors ceig)
                       (c/diag (:ivalues ceig))
                       (c/t (:vectors ceig))))))

;; SVD
;;
;; Sometimes this doesn't hold when Q has negative entries. Why not?
;; Jblas problem?  Or is my expectation wrong?
(let [Z (c/rand 4 7)
      Q (c/rnorm 4 7)]
  (let [{left :left values :values right :right} (c/svd Z)]
    (expect Z (c/* left
                   (c/diag values)
                   (c/t right))))
  (let [{left :left values :values right :right} (c/svd Q)]
    (expect Q (c/* left
                   (c/diag values)
                   (c/t right)))))

;; TODO test full SVD

;; matrix powers
(expect (c/* B B) (c/pow B 2))

;; properties of special random matrices
(let [H (c/rreflection n)
      P (c/rspectral n)
      G (c/cholesky P)
      I (c/rspectral (repeat n (double 1)))]
  (expect (c/id n) (c/* H (c/t H)))
  (expect (partial every? pos?) (:values (c/eigen P)))
  (expect P (c/* (c/t G) G))
  (expect (c/id n) I))

;; Matrix Functions
(let [v (c/column (range 9))]

  (expect (Math/exp 7.0)    (c/get (c/exp v) 7 0))
  (expect (Math/abs 7.0)    (c/get (c/abs v) 7 0))
  (expect (Math/acos 1.0)   (c/get (c/acos v) 1 0))
  (expect (Math/asin 1.0)   (c/get (c/asin v) 1 0))
  (expect (Math/atan 7.0)   (c/get (c/atan v) 7 0))
  (expect (Math/cbrt 7.0)   (c/get (c/cbrt v) 7 0))
  (expect (Math/ceil 7.0)   (c/get (c/ceil v) 7 0))
  (expect (Math/cos 7.0)    (c/get (c/cos v) 7 0))
  (expect (Math/cosh 7.0)   (c/get (c/cosh v) 7 0))
  (expect (Math/exp 7.0)    (c/get (c/exp v) 7 0))
  (expect (Math/floor 7.0)  (c/get (c/floor v) 7 0))
  (expect (Math/log 7.0)    (c/get (c/log v) 7 0))
  (expect (Math/log10 7.0)  (c/get (c/log10 v) 7 0))
  (expect (Math/signum 7.0) (c/get (c/signum v) 7 0))
  (expect (Math/sin 7.0)    (c/get (c/sin v) 7 0))
  (expect (Math/sinh 7.0)   (c/get (c/sinh v) 7 0))
  (expect (Math/sqrt 7.0)   (c/get (c/sqrt v) 7 0))
  (expect (Math/tanh 7.0)   (c/get (c/tanh v) 7 0))

  ;; These are the mutating versions, so pass in a fresh
  ;; matrix each time
  (expect (Math/exp 7.0)    (c/get (c/exp! (c/column (range 9))) 7 0))
  (expect (Math/abs 7.0)    (c/get (c/abs! (c/column (range 9))) 7 0))
  (expect (Math/acos 1.0)   (c/get (c/acos! (c/column (range 9))) 1 0))
  (expect (Math/asin 1.0)   (c/get (c/asin! (c/column (range 9))) 1 0))
  (expect (Math/atan 7.0)   (c/get (c/atan! (c/column (range 9))) 7 0))
  (expect (Math/cbrt 7.0)   (c/get (c/cbrt! (c/column (range 9))) 7 0))
  (expect (Math/ceil 7.0)   (c/get (c/ceil! (c/column (range 9))) 7 0))
  (expect (Math/cos 7.0)    (c/get (c/cos! (c/column (range 9))) 7 0))
  (expect (Math/cosh 7.0)   (c/get (c/cosh! (c/column (range 9))) 7 0))
  (expect (Math/exp 7.0)    (c/get (c/exp! (c/column (range 9))) 7 0))
  (expect (Math/floor 7.0)  (c/get (c/floor! (c/column (range 9))) 7 0))
  (expect (Math/log 7.0)    (c/get (c/log! (c/column (range 9))) 7 0))
  (expect (Math/log10 7.0)  (c/get (c/log10! (c/column (range 9))) 7 0))
  (expect (Math/signum 7.0) (c/get (c/signum! (c/column (range 9))) 7 0))
  (expect (Math/sin 7.0)    (c/get (c/sin! (c/column (range 9))) 7 0))
  (expect (Math/sinh 7.0)   (c/get (c/sinh! (c/column (range 9))) 7 0))
  (expect (Math/sqrt 7.0)   (c/get (c/sqrt! (c/column (range 9))) 7 0))
  (expect (Math/tanh 7.0)   (c/get (c/tanh! (c/column (range 9))) 7 0)))

(expect 2.0498889512140543
        (c/trace (c/with-seed 42 (c/rand 5 5))))
