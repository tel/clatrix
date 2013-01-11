(ns clatrix.matrix-api-test
  (:use clojure.test)
  (:require [clatrix.core :as c]))

(def M (c/matrix [[1 2 3] [4 5 6] [7 8 9]]))
(def V (c/matrix [[1 2 3 4 5]]))

(deftest implementation-key-test
  (is (= :clatrix (implementation-key M))))

(deftest dimension-info-test
  (is (= 2 (dimensionality M)))
  (is (= 1 (dimensionality V))))
