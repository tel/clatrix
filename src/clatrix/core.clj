(ns clatrix.core
  (:refer-clojure :exclude [get set rand vector?])
  (:use [slingshot.slingshot :only [throw+]])
  (:import [org.jblas DoubleMatrix Eigen Solve Singular MatrixFunctions]
           [java.io Writer]))

;;; A wrapper object over DoubleMatrix
;;;
(defrecord Matrix [^DoubleMatrix me]
  Object
  (toString [^Matrix mat]
    (str (list 'matrix
               (vec (map vec (vec (.toArray2 (:me mat)))))))))

;;; Java interop
;;; 
(defmacro promote-cfun* [defname name fname]
  (let [n (gensym)
        m (gensym)]
    `(~defname ~name
       ([^long ~n] (Matrix. (~fname ~n)))
       ([^long ~n ~m] (Matrix. (~fname ~n ~m))))))

(defmacro promote-mfun* [defname name fname]
  (let [m (gensym)]
    `(~defname ~name [^Matrix ~m] (~fname (:me ~m)))))

;;; Bootstrap the basics
;;; 
(defn matrix? [m]
  (isa? (class m) Matrix))
(defn get [^Matrix m ^long r ^long c]
  (.get (:me m) r c))
(defn set [^Matrix m ^long r ^long c ^double e]
  (.put (:me m) r c e))

(promote-mfun* defn- ncols .columns)
(promote-mfun* defn- nrows .rows)
(defn size    [^Matrix m] [(nrows m) (ncols m)])
(defn vector? [^Matrix m] (some #(== 1 %) (size m)))
(defn row?    [^Matrix m] (== 1 (first (size m))))
(defn column? [^Matrix m] (== 1 (second (size m))))
(defn square? [^Matrix m] (reduce == (size m)))

;;; Matrix creation
;;;

;;; Direct imports
(defn column [^doubles seq]
  (Matrix. (DoubleMatrix. (into-array Double/TYPE seq))))
(defn matrix [seq-of-seqs]
  (Matrix. (DoubleMatrix. (into-array (map #(into-array Double/TYPE %) seq-of-seqs)))))

;;; Create a constant matrix of a particular size.
(defn constant
  ([^long n ^double c] 
     (Matrix.
      (doto (DoubleMatrix/ones n)
        (.muli c))))
  ([^long n ^long m ^double c] 
     (Matrix.
      (doto (DoubleMatrix/ones n m) 
        (.muli c)))))

;;; Specific kinds of constant matrices and vectors
(promote-cfun* defn  ones DoubleMatrix/ones)
(promote-cfun* defn- zeros DoubleMatrix/zeros)

(defn diag [seq-or-matrix]
  (if (matrix? seq-or-matrix)
    (let [mat ^Matrix seq-or-matrix]
      ;; We'll extract largest diagonals from non-square matrices
      ;; since this isn't really a matrix algebraic property
      (let [n (apply min (size mat))]
        (map #(get mat % %) (range n))))
    (let [di ^doubles (seq seq-or-matrix)]
      (Matrix. (DoubleMatrix/diag (column di))))))

(defn id [^long n] (Matrix. (DoubleMatrix/eye n)))

;;; Matrix transformations
(defn t [^Matrix mat] (Matrix. (.transpose (:me mat))))

;;; Random matrices
(promote-cfun* defn  rand   DoubleMatrix/rand)
(promote-cfun* defn- randn* DoubleMatrix/randn) 

(defn rnorm
  ([^double mu ^double sigma ^long n ^long m]
     (Matrix.
      (doto (:me (randn* n m))
        (.muli sigma)
        (.addi mu))))
  ([^double mu ^double sigma ^long n]
     (Matrix.
      (doto (:me (randn* n))
        (.muli sigma)
        (.addi mu))))
  ([^long n ^long m] (randn* n m))
  ([^long n] (randn* n)))

;;; Stacking
;;;
(defn hstack [& vec-seq]
  (let [row-counts (map #(first (size %)) vec-seq)
        rows (first row-counts)]
    (if (every? (partial == rows) row-counts)
      (Matrix. (reduce #(DoubleMatrix/concatHorizontally %1 (:me %2))
                       (:me (first vec-seq))
                       (rest vec-seq))))))

(defn vstack [& vec-seq]
  (let [col-counts (map #(second (size %)) vec-seq)
        cols (first col-counts)]
    (if (every? (partial == cols) col-counts)
      (Matrix. (reduce #(DoubleMatrix/concatVertically %1 (:me %2))
                       (:me (first vec-seq))
                       (rest vec-seq))))))

;;; Matrix queries
;;;
(defn rows
  ([^Matrix mat]
     (let [[n m] (size mat)]
       (rows mat (range n))))
  ([^Matrix m ^longs idxs]
     (map #(Matrix. (.getRow (:me m) %)) idxs)))

(defn cols
  ([^Matrix mat]
     (let [[n m] (size mat)]
       (cols mat (range m))))
  ([^Matrix m ^longs idxs]
     (map #(Matrix. (.getColumn (:me m) %)) idxs)))

(defn as-vec
  "Converts a matrix object into a lazy sequence of its elements, row-major."
  [^Matrix m]
  (if (vector? m)
    (vec (.toArray (:me m)))
    (vec (map vec (vec (.toArray2 (:me m)))))))

;;; Slicing
;;;
(defn- iswild [sym] (= (str sym) "_"))

(defn- slicer
  ([^Matrix matrix rowspec colspec]
     (cond (and (iswild rowspec) (iswild colspec)) matrix
           (iswild rowspec) `(Matrix. (.getColumn (:me ~matrix) ~colspec))
           (iswild colspec) `(Matrix. (.getRow    (:me ~matrix) ~rowspec))
           :else            `(get                 ~matrix ~rowspec ~colspec)))
  ([^DoubleMatrix matrix rowspec colspec values]
     (let [m (gensym)
           form (cond (and (iswild rowspec) (iswild colspec)) `(.copy ~m ~values)
                      (iswild rowspec) `(.putColumn (:me ~m) ~colspec ~values)
                      (iswild colspec) `(.putRow    (:me ~m) ~rowspec ~values)
                      :else            `(set        ~m ~rowspec ~colspec ~values))]
       `(let [~m ~matrix]
          (do ~form ~m)))))

(defmacro slice
  "Slice is the primary function for accessing and modifying a matrix
at the single row, column, entry, or full matrix level. The
row/colspec variables are either an integer or the atom '_ signifying
that the index should run over all possible values for the row or
column index. If a fourth argument is passed it is assumed to be a
size-conforming entry, row, or matrix to be inserted into the spec'd
location."
  [^Matrix matrix rowspec colspec & values?]
  (apply slicer matrix rowspec colspec values?))

(defmacro slices
  "Identical interface to `slice` except that it returns a sequence"
  [^Matrix matrix rowspec colspec & values?]
  `(as-vec ~(apply slicer matrix rowspec colspec values?)))

(defn trace
  [^Matrix mat]
  (if (square? mat)
    (let [[n _] (size mat)]
      (reduce #(+ (slice mat %2 %2) %1) 0 (range n)))
    (throw+ {:error "Cannot take trace of non-square matrix."})))

;;; Hinting
;;;
(defn symmetric
  "Asserts that a matrix is symmetric."
  [^Matrix m] (with-meta m {:symmetric true}))

(defn positive
  "Asserts that a matrix is positive definite."
  [^Matrix m] (with-meta (symmetric m) {:positive true}))

(defn arbitrary
  "Asserts that a matrix is just arbitrary."  [^Matrix m] (with-meta m
  {:symmetric false :positive false}))

(defn symmetric? [^Matrix m] (:symmetric (meta m)))
(defn positive?  [^Matrix m] (:positive (meta m)))
(defn arbitrary? [^Matrix m] (not (or (symmetric? m) (positive? m))))

(defn maybe-symmetric
  "Checks to see if a matrix is symmetric, then hints if true."
  [^Matrix m]
  (if (or (symmetric? m) (= (t m) m))
    (symmetric m)
    (with-meta m {:symmetric false})))

(defn solve
  "Solves the equation AX = B for X. Flags include :pd :positive
  and :symmetric ensuring that optimized methods can be used"
  [^Matrix A ^Matrix B]
  (cond
   (psd? A)       (Solve/solvePositive (:me A) (:me B))
   (symmetric? A) (Solve/solveSymmetric (:me A) (:me B))
   :else          (Solve/solve (:me A) (:me B))))

(defn eigen
  "Computes the eigensystem (or generalized eigensystem) for a square
  matrix A. (eigen A :symmetric) uses optimized routines for symmetric
  matrices while (eigen A :symmetric B) computes the generalized eigenvectors x
  such that A x = L B x for symmetric A and B."
  ([^DoubleMatrix A]
     (cond
      (symmetric? A) (let [[vecs vals] (map #(Matrix. %)
                                            (seq (Eigen/symmetricEigenvectors (:me A))))]
                       {:vectors (cols vecs) :values (diag vals)})
      :else          (let [[vecs vals] (seq (Eigen/eigenvectors (:me A)))
                           rvecs (Matrix. (.real vecs))
                           ivecs (Matrix. (.imag vecs))
                           rvals (diag (Matrix. (.real vals)))
                           ivals (diag (Matrix. (.imag vals)))
                           out {:vectors (cols rvecs)
                                :values  rvals}]
                       (if (some (partial not= 0.0) ivals)
                         (merge out
                                {:ivectors (cols ivecs)
                                 :ivalues  ivals})
                         out))))
  ([^DoubleMatrix A ^DoubleMatrix B]
     (let [A (maybe-symmetric A)
           B (maybe-symmetric B)]
       (if (and (symmetric? A) (symmetric? B))
         (let [[vecs vals]
               (map #(Matrix. %)
                    (seq (Eigen/symmetricGeneralizedEigenvectors (:me A) (:me B))))]
           {:vectors (cols vecs) :values (as-vec vals)})
         (throw+ {:error "Cannot do generalized eigensystem for non-symmetric matrices."})))))

(defn maybe-positive
  "Checks to see if a matrix is symmetric, then hints if true."
  [^Matrix m]
  (let [m (check-symmetric m)]
    (if (symmetric? m)
      (let [{vals :values} (eigen m)]
        (if (every? #(> 0 %) vals)
          (positive m)
          m))
      m)))

;;; Viewing
;;;

(defmethod print-method Matrix [mat ^Writer w]
  (.write w "#<Matrix ")
  (.write w (str (size mat)))
  (.write w ">"))

;;; This one is pretty ugly, but whatever.
(defn pp
  "Pretty printer for matrices. The second and third optional
arguments are for large, precise matrices, specifying the amount of
elements to show and their precision respectively."
  ([^Matrix mat] (pp mat 3 4))
  ([^Matrix mat prec] (pp mat 3 prec))
  ([^Matrix mat nbits prec]
     (let [[n m] (size mat)
           small-rows (< n (* nbits 2))
           small-cols (< m (* nbits 2))
           rowset (if (< n (* nbits 2))
                    (range n)
                    (concat (range nbits) (range (- n nbits) n)))
           colset (if (< m (* nbits 2))
                    (range m)
                    (concat (range nbits) (range (- m nbits) m)))
           submat (apply hstack (cols (apply vstack (rows mat rowset)) colset))
           [n m] (size submat)
           fmt (str "% ." prec "e ")
           header (apply format " A %dx%d matrix" (size mat))]
       ;; Print the header
       (println header)
       (print " ")
       (doseq [i (range (dec (count header)))] (print "-"))
       (print "\n")
       ;; Print the matrix
       (if small-rows
         (doseq [i (range n)]
           (if small-cols
             (doseq [j (range m)]
               (print (format fmt (slice submat i j))))
             (do (doseq [j (range nbits)]
                   (print (format fmt (slice submat i j))))
                 (print " . ")
                 (doseq [j (range nbits (* 2 nbits))]
                   (print (format fmt (slice submat i j))))))
           (print "\n"))
         (do (doseq [i (range nbits)]
               (if small-cols
                 (doseq [j (range m)]
                   (print (format fmt (slice submat i j))))
                 (do (doseq [j (range nbits)]
                       (print (format fmt (slice submat i j))))
                     (print " . ")
                     (doseq [j (range nbits (* 2 nbits))]
                       (print (format fmt (slice submat i j))))))
               (print "\n"))
             (println " ... ")
             (doseq [i (range nbits (* 2 nbits))]
               (if small-cols
                 (doseq [j (range m)]
                   (print (format fmt (slice submat i j))))
                 (do (doseq [j (range nbits)]
                       (print (format fmt (slice submat i j))))
                     (print " . ")
                     (doseq [j (range nbits (* 2 nbits))]
                       (print (format fmt (slice submat i j))))))
               (print "\n")))))))