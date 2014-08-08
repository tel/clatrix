(ns clatrix.loading-test
  (:use [clojure.core.matrix]))

(defn foo []
  (doall (map deref (for [i (range 10)] (future (matrix :clatrix [0 1]))))))