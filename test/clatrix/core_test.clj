(ns clatrix.core-test
  (:use clojure.test)
  (:require [clatrix.core :as c]))

;;; # Utilities
;;;
(defmacro with-randmat [name & body]
  `(let [~name (c/rand (Math/floor (+ 10 (rand 100)))
                       (Math/floor (+ 10 (rand 100))))]
     ~@body))

;;; # Matrix creation
;;;

(deftest dense-matrix-inversion
  (with-randmat A
    (is (= A (c/matrix (c/dense A))))))

(deftest diag-involution
  (with-randmat A
    (is (= (c/diag A) (c/diag (c/diag (c/diag A)))))))

;;; ## Random matrices
;;;

;;; ## Element algebra
;;;

(deftest t-involutive
  (with-randmat A
    (is (= A (c/t (c/t A))))))

(deftest stack-explode-inversion
  (with-randmat A
    (is (= A (apply c/hstack (c/cols A))))
    (is (= A (apply c/vstack (c/rows A))))))

(deftest permute-reverse-involution
  (with-randmat A
    (let [[n m] (c/size A)
          rowspec (reverse (range n))
          colspec (reverse (range m))]
      (is (= A
             (c/permute
              (c/permute A
                         :rowspec rowspec
                         :colspec colspec)
              :rowspec rowspec
              :colspec colspec))))))

(deftest block-identity
  (let [A (c/id 30)
        B (c/id 10)]
    (is (= (c/id 30)
           (c/block [[B . .]
                     [. B .]
                     [. . B]])))))

(deftest block-reassembly
  (with-randmat A
    (let [[c1 c2 c3] (c/cols A [0 1 2])
          Ac (c/hstack c1 c2 c3)
          [r1 r2 r3] (c/rows A [0 1 2])
          Ar (c/vstack r1 r2 r3)]
      (is (= Ac (c/block [[c1 c2 c3]])))
      (is (= Ar (c/block [[r1]
                          [r2]
                          [r3]]))))))