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
    (is (c/clatrix? m))
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

(deftest maths-tests
  (let [m (m/matrix :clatrix [[1 2] [3 4]])
        v (m/matrix :clatrix [10 20])]
    (is (m/numerical? m))
    (is (m/equals m (m/add m 0)))
    (is (m/equals m (m/mul m 1)))
    (is (m/equals [50 110] (m/mmul m v)))
    (is (m/equals [[7.0 10.0] [15.0 22.0]] (m/mmul m m)))
    (is (m/equals 500 (m/mmul v v)))
    (is (m/equals [[1.0 (/ 2.0)] [(/ 3.0) (/ 4.0)]] (m/div m)))
    (is (m/equals [[1 1] [1 1]] (m/div m m)))))

(deftest map-tests
  (let [v (c/vector [1 2 3])
        v2 (c/vector [10 20 30])]
    (is (= (c/size v) (c/size v2) [3]))
    (is (m/equals [2 4 6] (m/emap clojure.core/+ v v)))
    (m/emap! clojure.core/+ v v2)
    (is (m/equals [11 22 33] v))))

(deftest matrix-tests
  (let [m (m/matrix :clatrix [[1 2 3] [3 4 5]])]
    (is (== 2 (m/dimensionality m)))
    (is (== 2 (m/dimension-count m 0)))
    (is (== 3 (m/dimension-count m 1)))
    (is (m/equals [1 3] (first (m/columns m))))
    (is (m/equals [1 2 3] (first (m/rows m))))))

(deftest construction-tests
  (let [m (c/matrix [[1 2] [3 4]])
        v (c/vector [1 2 3])]
    (is (= [2 2] (m/shape m)))
    (is (= [3] (m/shape v)))
    (is (== 10.0 (m/esum m)))
    (is (== 6.0 (m/esum v)))))

(deftest sequence-tests
  (is (= [1.0 2.0 3.0 4.0] (m/eseq (c/matrix [[1 2] [3 4]])))))

(deftest instance-tests
  (comp/instance-test (c/matrix [[1 2] [3 4]]))
  (comp/instance-test (c/matrix [[1 2] [3 4] [5 6]]))
  (comp/instance-test (c/vector [1 2]))
)

(deftest compliance-test
  (comp/compliance-test (c/matrix [[1 2] [3 4]]))
)

(deftest norm-tests
  (let [m (m/matrix :clatrix [[1 2] [3 4]])
        v (m/matrix :clatrix [1 2])]
    (is (m/equals (p/norm v 1) 3.0))
    (is (m/equals (p/norm v 2) 2.23606797749979))
    (is (m/equals (p/norm v 3) 2.080083823051904))
    (is (m/equals (p/norm v Double/POSITIVE_INFINITY) 2.0))
    (is (m/equals (p/norm m 1) 6.0))
    (is (m/equals (p/norm m 2) 5.464985704219043))
    (is (m/equals (p/norm m Double/POSITIVE_INFINITY) 7.0))))

(deftest decompositions-tests
  (let [m (m/matrix :clatrix [[2 -1 0] [-1 2 -1] [0 -1 2]])]
    (testing "QR decomposition"
        (let [{:keys [Q R]} (p/qr m {:return [:Q :R]})]
          (is (m/equals Q (m/matrix :clatrix
                                    [[-0.89442719 -0.35856858 0.26726124]
                                     [0.4472136 -0.71713717 0.53452248]
                                     [0 0.5976143 0.80178373]])))
          (is (m/equals R (m/matrix :clatrix
                                    [[-2.23606798 1.78885438 -0.4472136]
                                     [ 0 -1.67332005 1.91236577]
                                     [ 0 0 1.06904497]])))
          (is (m/equals (m/mmul Q R) m))))
    (testing "Cholesky decomposition"
      (let [{:keys [L L*]} (p/cholesky m {:return [:L :L*]})]
        (is (m/equals L (m/matrix :clatrix
                                  [[1.41421356 0 0]
                                    [-0.70710678 1.22474487 0]
                                    [0 -0.81649658 1.15470054]])
                      1e-5))
        (is (m/equals L* (m/matrix :clatrix
                                  [[1.41421356 -0.70710678 0]
                                   [0 1.22474487 -0.81649658]
                                   [0 0 1.15470054]])
                      1e-5))))

    (testing "LU decomposition"
      (let [{:keys [P L U]} (p/lu m {:return [:P :L :U]})]
        (is (m/equals P (m/identity-matrix :clatrix 3)))
        (is (m/equals L (m/matrix :clatrix
                                  [[1 0 0]
                                   [-0.5 1 0]
                                   [0 -0.66667 1]])))
        (is (m/equals U (m/matrix :clatrix
                                  [[2 -1 0]
                                   [0 1.5 -1]
                                   [0 0 1.33333]])))))

    (testing "SVD decomposition"
      (let [{:keys [U S V*]} (p/svd m {:return [:U :S :V*]})]
        (is (m/equals U (m/matrix :clatrix
                                  [[-0.5 -0.707106 0.5]
                                   [0.707106 0 0.707106]
                                   [-0.5 0.707106 0.5]])
                      1e-5))
        (is (m/equals S (m/matrix :clatrix [3.4142135 2 0.5857864])
                      1e-5))
        (is (m/equals V* (m/matrix :clatrix
                                   [[-0.5 -0.707106 0.5]
                                    [0.707106 0 0.707106]
                                    [-0.5 0.707106 0.5]])
                      1e-5))))

    (testing "Eigen decomposition"
      (let [{:keys [Q rA iA]} (p/eigen m {:return [:Q :rA :iA]})]
        (is (m/equals Q (m/matrix :clatrix
                                  [[-0.5 -0.707106 0.5]
                                   [0.707106 0 0.707106]
                                   [-0.5 0.707106 0.5]]) 1e-5))
        (is (m/equals rA (m/matrix :clatrix [3.41421 0 0 0 2 0 0 0 0.58578])
                      1e-5))
        (is (m/equals iA (m/new-vector :clatrix 9)))))

    (testing "Symmetric Eigen decomposition"
      (let [sym (m/matrix :clatrix [[1 3 0]
                                    [3 2 6]
                                    [0 6 5]])
            {:keys [Q rA iA]} (p/eigen sym {:return [:Q :rA :iA]
                                            :symmetric true})]
        (is (m/equals Q (m/matrix :clatrix
                                  [[-0.45327 0.86657 0.20878]
                                   [0.73883 0.23422 0.63186]
                                   [-0.49865 -0.44066 0.74642]])
                      1e-5))
        (is (m/equals rA (m/matrix :clatrix [-3.88999 0 0 0 1.81086 0 0 0 10.07912])
                      1e-5))
        (is (= iA nil))))))

(deftest solve-tests
  (let [x (m/matrix :clatrix [[3 2 -1]
                              [2 -2 4]
                              [-1 0.5 -1]])
        y (m/matrix :clatrix [1 -2 0])]
    (is (m/equals (m/matrix :clatrix [[1] [-2] [-2]])
                  (p/solve x y)))))

(deftest least-squares-tests
  (let [x (m/matrix :clatrix [[3 2 -1]
                              [2 -2 4]
                              [-1 0.5 -1]])
        y (m/matrix :clatrix [1 -2 0])]
    (is (m/equals (m/matrix :clatrix [[1] [-2] [-2]])
                  (p/least-squares x y)))))

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
