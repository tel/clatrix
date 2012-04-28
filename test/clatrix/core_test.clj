(ns clatrix.core-test
  (:use expectations)
  (:require [clatrix.core :as c]))

(let [A (c/rand 10 15)
      S (c/rand 10 10)
      [n m] (c/size A)
      [ns _] (c/size S)
      ridx (range n)
      cidx (range m)]
  ;; creation
  (expect A (c/matrix (c/dense A)))

  ;; structure/extraction
  (expect (c/diag A) (c/diag (c/diag (c/diag A))))

  ;; algebra
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
                         [r3]]))))