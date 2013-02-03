(ns clatrix.matrix-api-test
  (:use clojure.test
        core.matrix.protocols)
  (:require [clatrix.core :as c]
            [core.matrix :as m]
            [core.matrix.protocols :as p]
            [core.matrix.compliance-tester :as comp]))

(deftest compliance-test
  (comp/compliance-test (c/matrix [[1 2] [3 4]])))
