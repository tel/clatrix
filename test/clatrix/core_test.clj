(ns clatrix.core-test
  (:use expectations)
  (:import [clatrix.core Matrix]
           [java.io StringReader PushbackReader])
  (:require [clatrix.core :as c]))

(defn read* [str]
  (read (PushbackReader. (StringReader. str))))

(defn truthy? [v]
  (not (or (false? v) (nil? v))))

(let [[n m] [10 15]
      ns n
      A (c/rnorm n m)
      B (c/* (c/t A) A)
      S (c/rnorm ns ns)
      p (c/rnorm m)
      q (c/rnorm ns)
      ridx (range n)
      cidx (range m)]

  ;; properties of A/S
  (expect Matrix A)
  (expect Matrix S)
  (expect Matrix p)
  (expect Matrix q)
  
  (given A
         (expect c/size [n m]
                 c/matrix? true
                 c/square? false
                 c/column? false
                 c/vector? false
                 c/row?    false))
  (given S
         (expect c/size [ns ns]
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
         (expect c/size [ns 1]
                 c/matrix? true
                 c/square? false
                 c/column? true
                 c/vector? true
                 c/row?    false))

  (let [z (rand)]
    (c/set A 0 0 z)
    (expect z (c/get A 0 0)))

  (expect `(c/matrix ~(map #(map double %) [[1 2] [3 4]]))
          (read* (str (c/matrix [[1 2] [3 4]]))))
  
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
  (expect A (c/hstack (c/cols A)))
  (expect A (c/vstack (c/rows A)))

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
    ;;   (V + iV I)(L + iL I)(Vt - iVt I)
    ;; = I (iV iL iVt + iV  L  Vt +  V iL  Vt -  V  L iVt)
    ;; +   (V   L  Vt +  V iL iVt + iV  L iVt - iV iL  Vt)
    ;;
    ;; where `iA` is the imaginary part of some matrix with real part
    ;; `A`.
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
  ;; Doesn't hold if Z has negative entries. Why not?
  (let [Z (c/rand 4 7)]
    (let [{left :left values :values right :right} (c/svd Z)]
      (expect #(> 1.0E-12 (Math/abs %))
              (- (c/norm Z)
                 (c/norm (c/* left
                              (c/diag values)
                              (c/t right)))))))
  
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
    (expect (c/id n) I)))