(ns clatrix.core
  (:refer-clojure :exclude [get set rand vector?])
  (:use [slingshot.slingshot :only [throw+]])
  (:import [org.jblas DoubleMatrix Eigen Solve Singular MatrixFunctions]
           [java.io Writer]))

;;; Clatrix is a fast matrix library for Clojure written atop JBlas'
;;; ATLAS/LAPACK bindings. It's not intended to be the alpha and omega
;;; of linear algebra hosted on Clojure, but it should provide most of
;;; the basics for enabling more scientific/mathematical computation.

;;; # The clatrix matrix

;;; Matrices are implemented as a thin wrapper over JBlas'
;;; `DoubleMatrix` class. The wrapper is useful to program metadata
;;; and protocols into the Matrix implementation and also to hide the
;;; underlying Java methods. It's not hard to access them (they're
;;; available through the `:me` keyword), but their use is clearly
;;; dissuaded.
(defrecord Matrix [^DoubleMatrix me]
  Object
  (toString [^Matrix mat]
    (str (list 'matrix
               (vec (map vec (vec (.toArray2 (:me mat)))))))))

;;; # Java interop
;;; 
;;; Clatrix lifts a lot of methods directly out of JBlas. Here are a
;;; few convenience macros used for the thinnest parts of the wrapper.
(defmacro promote-cfun* [defname name fname]
  (let [n (gensym)
        m (gensym)]
    `(~defname ~name
       ([^long ~n] (Matrix. (~fname ~n)))
       ([^long ~n ~m] (Matrix. (~fname ~n ~m))))))

(defmacro promote-mfun* [defname name fname]
  (let [m (gensym)]
    `(~defname ~name [^Matrix ~m] (~fname (:me ~m)))))

;;; # Basics of matrix objects
;;;
;;; In linear algebra, matrices are two-dimensional arrays of
;;; doubles. The object `Matrix` is our particular instantiation.

(defn matrix? [m]
  (isa? (class m) Matrix))

;;; The most fundamental question about a matrix is its size. This
;;; also defines a number of other ideas such as whether a matrix is a
;;; column or row vector. Columns are default, though, by convention,
;;; matrices are sometimes represented as nested seqs in row-major
;;; order.

(promote-mfun* defn- ncols .columns)
(promote-mfun* defn- nrows .rows)
(defn size    [^Matrix m] [(nrows m) (ncols m)])
(defn vector? [^Matrix m] (some #(== 1 %) (size m)))
(defn row?    [^Matrix m] (== 1 (first (size m))))
(defn column? [^Matrix m] (== 1 (second (size m))))
(defn square? [^Matrix m] (reduce == (size m)))

;;; The most basic matrix operation is elementwise getting and
;;; setting; setting should be dissuaded as well for a Clojure
;;; wrapper, but it's too useful to hide.

(defn get [^Matrix m ^long r ^long c]
  (.get (:me m) r c))
(defn set [^Matrix m ^long r ^long c ^double e]
  (.put (:me m) r c e))

;;; Already this is sufficient to get some algebraic matrix properties

(defn trace
  "Computes the trace of a matrix, the sum of its diagonal elements."
  [^Matrix mat]
  (if (square? mat)
    (let [[n _] (size mat)]
      (reduce #(+ (get mat %2 %2) %1) 0 (range n)))
    (throw+ {:error "Cannot take trace of non-square matrix."})))

;;; We can also map the entire matrices back into Clojure data
;;; structures like 2-nested vectors.

(defn dense
  "Converts a matrix object into a lazy sequence of its element,
  row-major."
  [^Matrix m]
  (vec (map vec (vec (.toArray2 (:me m))))))

(defn as-vec
  "Converts a matrix object into a lazy sequence of its elements,
  row-major. Treats `vector?` type matrices specially, though, and
  flattening the return to a single vector."
  [^Matrix m]
  (if (vector? m)
    (vec (.toArray (:me m)))
    (vec (map vec (vec (.toArray2 (:me m)))))))

;;; # Matrix creation
;;;
;;; Matrices can be created from a number of simple specifications,
;;; such as (1) direct coercions, (2) element-constant matrices and
;;; vectors, and (3) identity matrices.

(defn column
  "Creates a column Matrix from a seq of its elements."
  [^doubles seq]
  (Matrix. (DoubleMatrix. (into-array Double/TYPE seq))))

(defn matrix
  "Creates a Matrix from a seq of seqs, specifying the matrix in
  row-major order. The length of each seq must be identical."
  [seq-of-seqs]
  (let [lengths (map count seq-of-seqs)
        l0      (first lengths)]
    (if (every? (partial = l0) lengths)
      (Matrix. (DoubleMatrix. (into-array (map #(into-array Double/TYPE %) seq-of-seqs))))
      (throw+ {:error "Cannot create a ragged matrix."}))))

(defn diag
  "Creates a diagonal matrix from a seq or pulls the diagonal of a
  matrix out as a seq."
  [seq-or-matrix]
  (if (matrix? seq-or-matrix)
    (let [mat ^Matrix seq-or-matrix]
      ;; We'll extract largest diagonals from non-square matrices
      ;; since this isn't really a matrix algebraic property
      (let [n (apply min (size mat))]
        (map #(get mat % %) (range n))))
    (let [di ^doubles (seq seq-or-matrix)]
      (Matrix. (DoubleMatrix/diag (column di))))))

(defn constant
  "Create a column or matrix filled with a constant value."
  ([^long n ^double c] 
     (Matrix.
      (doto (DoubleMatrix/ones n)
        (.muli c))))
  ([^long n ^long m ^double c] 
     (Matrix.
      (doto (DoubleMatrix/ones n m) 
        (.muli c)))))

(promote-cfun* defn  ones DoubleMatrix/ones)
(promote-cfun* defn- zeros DoubleMatrix/zeros)

(defn id
  "The `n`x`n` identity matrix."
  [^long n] (Matrix. (DoubleMatrix/eye n)))

;;; ## Random matrices
;;;
;;; It's also useful to generate random matrices. These are
;;; elementwise independent Unif(0,1) and normal.

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

;;; ## Element algebra
;;;
;;; Matrices can also be permuted, flipped, transposed, stacked, and
;;; split to form more complex matrices. This "element algebra" is a
;;; powerful way of building more complex matrices.

(defn t
  "The transpose of a matrix."
  [^Matrix mat] (Matrix. (.transpose (:me mat))))

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

(defn rows
  "Breaks a matrix into its constituent rows. By default, all the rows
are returned in normal order, but indices can also be specified."
  ([^Matrix mat]
     (let [[n m] (size mat)]
       (rows mat (range n))))
  ([^Matrix m ^longs idxs]
     (map #(Matrix. (.getRow (:me m) %)) idxs)))

(defn cols
  "Breaks a matrix into its constituent columns. By default, all the
columns are returned in normal order, but indices can also be
specified."
  ([^Matrix mat]
     (let [[n m] (size mat)]
       (cols mat (range m))))
  ([^Matrix m ^longs idxs]
     (map #(Matrix. (.getColumn (:me m) %)) idxs)))

(defn permute
  "Permutes the rows and the columns of a matrix"
  [^Matrix mat & {:keys [rowspec colspec]}]
  (let [[n m] (size mat)]
    (cond (and rowspec (some #(> % (dec n)) rowspec))
          (throw+ {:error "Row index out of bounds" :num-rows n :rowspec rowspec})
          (and colspec (some #(> % (dec m)) colspec))
          (throw+ {:error "Column index out of bounds" :num-columns n :colspec rowspec})

          :else
          (do
            (let [mat1 (if rowspec
                         (apply vstack (rows mat rowspec))
                         mat)
                  mat2 (if colspec
                         (apply hstack (cols mat1 colspec))
                         mat1)]
              mat2)))))

;;; ### Block matrices
;;;
;;; Block matrix syntax is a very convenient way of building larger
;;; matrices from smaller ones. Clatrix implements block matrix syntax
;;; in a convenient manner.
;;;
;;;     (block [[A 1 1 0]
;;;             [0 B . .]
;;;             [_ _ C .]
;;;             [_ _ D E]])
;;;
;;; Clatrix uses size constraint propagation to determine the proper
;;; sizes of the constant and 0 matrices.

(defn- iswild
  "Defines a wildcard symbol, used in `block`, `slice`, and `slices`."
  [sym]
  ;; This is a sort of silly way to do it, but I can't get the regex
  ;; to work for both '_ and #'user/_
  (let [name1 (second (re-find #"/(.+)" (str sym)))
        name2 (str sym)]
    (or (= name1 "_")
        (= name1 ".")
        (= name1 "*")
        (= name2 "_")
        (= name2 ".")
        (= name2 "*"))))

(defn- make-constr
  "Subfunction for `block`. Makes a size-constraint hash for building
  block matrices."
  [e]
  (if (matrix? e)
    {:matrix e
     :rows (first (size e))
     :cols (second (size e))}
    {:constant e}))

(defn- update-hash-with-constr
  "Examines a new constraint against an old block-hash, `hsh`. Updates the
  hash at position [i j] to respect the current constraint `constr`
  according to `key`"
  [hsh constr i j key]
  (let [{n key} constr
        old (hsh [i j])]
    (if (key old)
      (if (not= (key old) n) ;the constraint doesn't match the hash, uh oh
        (throw+ {:error "Block matrix diagram sizes are inconsistent."
                 :type :constraint-error
                 :location [i j]})
        hsh)
      ;; if there isn't an old key then we can fix that constraint now
      (assoc hsh [i j]
             (assoc old key n)))))

(defn- block-fn
  "Creates a block matrix. Any number `n` represents the all-`n`
  matrix of an appropriate size to make the matrix."
  [matrices]
  ;; We must do size-constraint propagation along the rows and columns
  ;; of the block-diagram in order to (a) ensure that the input isn't
  ;; in error and (b) find the proper sizes for the constant matrices.
  (let [n       (count matrices)
        lengths (map count matrices)
        m       (first lengths)]
    (if (every? (partial == m) lengths)
      
      ;; Build the constraints map
      (let [constrs (map #(map make-constr %) matrices)
            indices (for [i (range n) j (range m)] [i j])
            ;; The initial hash map contains what we know before
            ;; constraint propagation.
            init-map (reduce (fn [hsh [i j]]
                               (assoc hsh [i j]
                                      (nth (nth constrs i) j)))
                             (hash-map) indices)
            ;; Walk over the index set and propagate all the constraints
            ;; over each row and column
            constraint-map
            (reduce
             (fn [hash [i j]]
               (let [constr (nth (nth constrs i) j)]
                 (if (or (:rows constr) (:cols constr))
                   ;; Look up and to the left for constraint violations,
                   ;; locking in constraints if they don't already exist
                   (reduce #(update-hash-with-constr %1 constr i %2 :rows)
                           (reduce #(update-hash-with-constr %1 constr %2 j :cols)
                                   hash (range n))
                           (range m))
                   hash))) init-map indices)]
        ;; Use the constraint map to build the final matrix
        (apply vstack
               (for [i (range n)]
                 (apply hstack
                        (for [j (range m)]
                          (let [constr (constraint-map [i j])]
                            (if (:matrix constr)
                              (:matrix constr)
                              ;; Constants are assumed to be 1x1
                              ;; unless otherwise constrained
                              (constant (:rows constr 1)
                                        (:cols constr 1)
                                        (:constant constr)))))))))
      (throw+ {:error "Block matrices cannot be ragged."}))))

(defmacro block
  "Creates a block matrix using normal block matrix syntax written as
  a row-major ordered vector of vectors. Each entry in the
  specification can be a `Matrix`, a number, or a null symbol (either
  `.` or `_`). Numbers are translated as constant matrices of the
  appropriate size to complete the block matrix. Null symbols are
  considered as constant 0 matrices and are also automatically
  constrained to be the proper size. Any integers which do not share a
  row or a column with a larger matrix are assumed to be 1x1 sized."
  [blockspec]
  `(block-fn ~(vec (map #(vec (map (fn [e] (if (iswild e) 0 e)) %)) blockspec))))


;;; # Slicing 
;;;
;;; As a more convenient API than looping over `get` and `set`, we
;;; have slice notation. This uses wildcard symbols (like in `block`)
;;; in order to represent full index sets.
;;;
;;;     (slice A _ 5) ; ==> the whole 5th column
;;;     (slice A 4 5) ; ==> (get A 4 5)
;;;
;;; The slice macro also overloads setters. For instance
;;;
;;;     (slice A _ 1 (column (range 10)))
;;;
;;; replaces the 2nd column of `A` with `[0 1 2 3 4 5 6 7 8 9]`.

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

;;; # Hinting
;;;
;;; Many more complex matrix operations can be specialized for certain
;;; kinds of matrices. In particular, symmetric and positive definite
;;; matrices are much easier to handle. Clatrix doesn't natrually know
;;; which matrices have these special properties, but hinting
;;; functions can be used to assert that certain matrices are indeed
;;; symmetric or positive definite.
;;;
;;; Most operations in Clatrix create new objects, thus the assertions
;;; do not propagate. If they are, however, they can be removed by
;;; asserting the matrix `arbitrary`.

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

;;; # Linear algebra
;;;

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

;;; # Viewing the matrices
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