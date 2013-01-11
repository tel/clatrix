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

(deftest indexed-access-test
  (is (thrown? UnsupportedOperationException (get-1d M 1)))
  (is (= 6.0 (get-2d M 1 2)))
  (is (thrown? UnsupportedOperationException (get-nd M [1 2 3]))))

(deftest indexed-setting-test
  (is (thrown? UnsupportedOperationException (set-1d M 2 9.0)))
  (is (= (c/matrix [[1 2 9 4 5]]) (set-2d V 0 2 9.0)))
  (is (thrown? UnsupportedOperationException (set-nd M [1 2 3] 9.0))))
