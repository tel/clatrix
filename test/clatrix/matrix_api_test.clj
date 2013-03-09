(ns clatrix.matrix-api-test
  (:use clojure.test
        clojure.core.matrix.protocols)
  (:require [clatrix.core :as c]
            [clojure.core.matrix :as m]
            [clojure.core.matrix.protocols :as p]
            [clojure.core.matrix.compliance-tester :as comp]))

(deftest compliance-test
  (comp/compliance-test (c/matrix [[1 2] [3 4]])))

(comment
  (def M (c/matrix [[1 2 3] [4 5 6] [7 8 9]]))
  (def v (c/vector [1 2 3]))
  (c/* v M)

  (c/vector? v)
  (c/column? v)
  (c/matrix? v)
  (p/is-scalar? (p/coerce 1))
  (c/get v 2)
  (p/get-1d v 2)
  (p/vector-dot v v)
  (def v (c/vector [1 2 3]))
  (p/length v)
  (c/t v)

  (p/length-squared v)
  (p/normalise v)

  (c/get (c/slice M 0 _) 2)
  (c/get (c/slice M _ 0) 2)

  (c/get M [0 4]))
