(ns clatrix.core-test
  (:use expectations)
  (:import [clatrix.core Matrix])
  (:require [clatrix.core :as c]))

(let [[n m] [10 15]
      ns n
      A (c/rand n m)
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
  (expect (double (* 7 20)) (c/trace (c/+ (c/id 20) (c/id 20) 5)))
  (expect A (c/- (c/+ A A) A)))