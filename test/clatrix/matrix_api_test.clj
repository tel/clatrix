(ns clatrix.matrix-api-test
  (:use clojure.test
        clojure.core.matrix.protocols)
  (:require [clatrix.core :as c]
            [clojure.core.matrix :as m]
            [clojure.core.matrix.protocols :as p]
            [clojure.core.matrix.compliance-tester :as comp]))


(deftest regressions
  (let [m (m/matrix :clatrix [1 2 3])
        wm (clojure.core.matrix.impl.wrappers/wrap-nd m)]
    (is (== 1 (dimensionality m)))
    (is (= [3] (seq (m/shape m))))
    (is (m/equals m wm))
    (is (= [3] (seq (m/shape wm))))
    (is (== 1 (m/dimensionality wm)))
    (is (p/is-vector? wm))
    (is (== 2.0 (m/mget wm 1)))
    (is (m/equals [2.0] (m/subvector wm 1 1)))
    (is (m/equals [2.0] (m/subvector m 1 1))))
  (is (m/equals [1.0 2.0] (m/subvector (m/matrix :clatrix [1 2 3]) 0 2)))
  (is (m/equals [1.0 2.0] (m/coerce [] (m/matrix :clatrix [1 2]))))
  (is (every? number? (m/slices (m/matrix :clatrix '(1 2 3)))))
  (is (= [1.0 2.0 3.0] (m/eseq (m/matrix :clatrix '(1 2 3)))))
  (is (= [1.0 2.0 3.0] (m/slices (m/matrix :clatrix [1 2 3]))))
  (is (= [1.0 2.0 3.0 4.0] (m/eseq (m/matrix :clatrix [[1 2] [3 4]]))))
  (is (m/equals [1 2 3] (m/matrix :clatrix [1 2 3]))))

(deftest matrix-tests
  (let [m (m/matrix :clatrix [[1 2 3] [3 4 5]])]
    (is (== 2 (m/dimensionality m)))
    (is (== 2 (m/dimension-count m 0)))
    (is (== 2 (m/dimension-count m 1)))
    (is (equals [1 3] (first (m/columns m))))
    (is (equals [1 2 3] (first (m/rows m))))))

(deftest compliance-test
  (comp/compliance-test (c/matrix [[1 2] [3 4]]))
)

(comment
  (def M (c/matrix [[1 2 3] [4 5 6] [7 8 9]]))  ;; 3x3 Matrix
  (get-shape M) ;; [3 3]
  (dimensionality M) ;; 2
  (get-shape M) ;; [3 3]

  (def r (get-row M 0)) ;; 1x3 Matrix
  (get-shape r)  ;; [1 3]
  (dimensionality r) ;; 2, as its a matrix.  Should this be one ?

  (get-column r 1) ;; 1x1 matrix: [2]
  (get-shape (get-column r 1)) ;; [1 1]
  (dimensionality (get-column r 1)) ;; 2.  As its still a matrix.  Should this be 0 ?
  (get-row r 0) ;; 1x3 matrix


  (cons M (c/matrix [[]]))
  (c/toString M)
  (count M)

  (element-seq r)

  (get-row r 0) ;; a 1x3 matrix
  (first (element-seq r))

  (def v (c/vector [1 2 3]))
  (def cv (c/matrix [1 2 3]))

  (some #(== 1 %) (c/size v))

  (c/normalize (c/matrix [3 4]))
  (element-seq cv)
  (clojure.core/flatten M)
  (tree-seq sequential? seq M)

  (c/size (c/matrix [2]))
  (== 1 1 1)

  (.vector? (c/t (c/matrix [1 2])))
  (.vector? (c/matrix (c/dotom .transpose (c/matrix [1 2]))))


  (m/vec?
   (clatrix.core/matrix [[1.0] [2.0]]))

  (dimensionality cv)

  (dimensionality cv)
  (dimensionality (first (get-major-slice-seq M)))
  (map #(c/slice M % _) (range (c/nrows M)))

  (m/slices M)
  (range (c/nrows M))

  (map get-0d (m/slices cv))

  (c/normalize (c/matrix [3 4]))

  (c/vector? cv)
  (c/set v 1 0 14)
  (c/set-1d (c/matrix M) 0 1 3)
  (m/dot v v)
  (length (c/matrix [3 4]))
  (equals [1 2 4] (m/mset (c/matrix [1 2 3]) 2 4))
  (c/matrix M)

  (filter (partial < 0) (c/size v))
  (m/slices M)
  (c/slices M _ _)
  (get-major-slice M 2)
  (get-major-slice-seq M)
  (dimensionality v)

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
  (c/get (c/slice M _ 0) 2))
