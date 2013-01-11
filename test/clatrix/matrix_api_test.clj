(ns clatrix.matrix-api-test
  (:use clojure.test
        core.matrix.protocols)
  (:require [clatrix.core :as c]))

(def M (c/matrix [[1 2 3] [4 5 6] [7 8 9] [10 11 12]]))
(def V (c/matrix [[1 2 3 4 5]]))
(def X (c/matrix [[1]]))

(deftest implementation-key-test
  (is (= :clatrix (implementation-key M))))

(deftest dimension-info-test
  (is (= 2 (dimensionality M)))
  (is (= 1 (dimensionality V)))
  (is (false? (is-scalar? M)))
  (is (true? (is-scalar? X)))
  (is (false? (is-vector? M)))
  (is (true? (is-vector? V)))
  (is (= 4 (dimension-count M 1)))
  (is (= 3 (dimension-count M 2)))) 
