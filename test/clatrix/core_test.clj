(ns clatrix.core-test
  (:use expectations)
  (:import [clatrix.core Matrix])
  (:require [clatrix.core :as c]))

(let [[n m] [10 15]
      ns n
      A (c/rand n m)
      B (c/* (c/t A) A)
      S (c/rand ns ns)
      ridx (range n)
      cidx (range m)]

  ;; properties of A/S
  (expect Matrix A)
  (expect Matrix S)
  (given A
         (expect c/size [n m]))
  (given S
         (expect c/size [ns ns]))

  ;; properties of id
  (given (c/id n)
         (expect c/size [n n]
                 c/trace (double n)))

  ;; conversion from Clojure types is invertible
  (expect A (c/matrix (c/dense A)))
  
  ;; diagonal structure becomes 2-parity involutive
  (expect (c/diag A) (c/diag (c/diag (c/diag A))))

  ;; structure algebraic constraints
  (expect A (c/t (c/t A)))
  (expect A (apply c/hstack (c/cols A)))
  (expect A (apply c/vstack (c/rows A)))

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

  ;; SVD decomposition
  (let [svd (c/svd A)]
    (expect A (c/* (:left svd)
                   (c/diag (:values svd))
                   (c/t (:right svd)))))
  
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