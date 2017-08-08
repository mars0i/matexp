(ns matexp.moorepenrose
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.linear :as lin]
            [clojure.core.matrix.utils :as u]))


;; Questions:
;; How should tolerance be calculated? 
;; Complex matrices are not fully part of core.matrix at present, right?
;; So I should use transpose rather than hermitian-transpose?
;; Eventually this needs to use hermitian-transpose on complex matrices.
(defn pinv 
  "Moore-Penrose pseudoinverse of matrix m calculated using svd.  tolerance 
  defaults to 0 during calculation. Absolute values of singular values that 
  are less than or equal to tolerance will be treated as zero.  Won't work
  if m's implementation doesn't have an implementation of svd unless svd
  is available for the current-implementation."
  ([m] (pinv m 0.0)) ;; default tolerance should be larger--maybe calculated
  ([m tolerance]
   (let [treat-as-nonzero? (fn [x] (> (Math/abs x) tolerance))
         diag-pinv (fn [singular-vals rows cols] ; pinv only for diagonal rectangular matrices
                     (let [smaller-dim (min rows cols)
                           sigma+ (m/ensure-mutable (m/zero-matrix cols rows))] ; dims transposed
                       (u/doseq-indexed [x singular-vals i]
                          (when (treat-as-nonzero? x)
                            (m/mset! sigma+ i i (/ x))))
                       sigma+))
         [rows cols] (m/shape m)
         {:keys [U S V*]} (lin/svd m)
         S+ (diag-pinv S rows cols)]
     (m/mmul (m/transpose V*) S+ (m/transpose U)))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; clatrix bug:
; 
; user=> (set-current-implementation :vectorz)
; :vectorz
; user=> (pm (pinv (matrix :vectorz [[1 2 3][4 5 6][7 8 9][10 11 12]]) 0.0000000000001))
; [[-0.483 -0.244 -0.006  0.233]
;  [-0.033 -0.011  0.011  0.033]
;  [ 0.417  0.222  0.028 -0.167]]
; nil
; user=> (pm (pinv (matrix :persistent-vector [[1 2 3][4 5 6][7 8 9][10 11 12]]) 0.0000000000001))
; [[-0.483 -0.244 -0.006  0.233]
;  [-0.033 -0.011  0.011  0.033]
;  [ 0.417  0.222  0.028 -0.167]]
; nil
; user=> (pm (pinv (matrix :ndarray [[1 2 3][4 5 6][7 8 9][10 11 12]]) 0.0000000000001))
; [[-0.483 -0.244 -0.006  0.233]
;  [-0.033 -0.011  0.011  0.033]
;  [ 0.417  0.222  0.028 -0.167]]
; nil
; user=> (pm (pinv (matrix :aljabr [[1 2 3][4 5 6][7 8 9][10 11 12]]) 0.0000000000001))
; [[-0.483 -0.244 -0.006  0.233]
;  [-0.033 -0.011  0.011  0.033]
;  [ 0.417  0.222  0.028 -0.167]]
; nil
; user=> (pm (pinv (matrix :clatrix [[1 2 3][4 5 6][7 8 9][10 11 12]]) 0.0000000000001))
; ExceptionInfo throw+: {:exception "Matrix products must have compatible sizes.", :a-cols 4, :b-rows 3}  slingshot.support/stack-trace (support.clj:201)
; 
; user=> (pst)
; ExceptionInfo throw+: {:exception "Matrix products must have compatible sizes.", :a-cols 4, :b-rows 3} {:exception "Matrix products must have compatible sizes.", :a-cols 4, :b-rows 3}
; 	slingshot.support/stack-trace (support.clj:201)
; 	clatrix.core/* (core.clj:1008)
; 	clatrix.core/* (core.clj:1002)
; 	clatrix.core/eval22125/fn--22220 (core.clj:1566)
; 	clojure.core.matrix.protocols/eval3411/fn--3425/G--3400--3432 (protocols.cljc:534)
; 	clojure.lang.ArraySeq.reduce (ArraySeq.java:109)
; 	clojure.core/reduce (core.clj:6544)
; 	clojure.core/reduce (core.clj:6527)
; 	clojure.core.matrix/mmul (matrix.cljc:1487)
; 	clojure.core.matrix/mmul (matrix.cljc:1473)
; 	matexp.moorepenrose/pinv (moorepenrose.clj:31)
; 	matexp.moorepenrose/pinv (moorepenrose.clj:13)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Verification of Moore-Penrose properties for one example:
; 
; user=> (def m (matrix :vectorz [[1 2 3][4 5 6][7 8 9][10 11 12]]))
; #'user/m
; user=> (def m+ (pinv (matrix :vectorz [[1 2 3][4 5 6][7 8 9][10 11 12]]) 0.000000000001))
; #'user/m+
; user=> (pm (sub m (mmul m m+ m)))
; [[0.000 0.000 -0.000]
;  [0.000 0.000 -0.000]
;  [0.000 0.000 -0.000]
;  [0.000 0.000 -0.000]]
; nil
; user=> (pm (sub m+ (mmul m+ m m+)))
; [[ 0.000  0.000  0.000 -0.000]
;  [ 0.000  0.000 -0.000 -0.000]
;  [-0.000 -0.000 -0.000  0.000]]
; nil
; user=> (let [c (mmul m m+)] (pm (sub c (transpose c))))
; [[ 0.000 0.000 -0.000 -0.000]
;  [-0.000 0.000 -0.000 -0.000]
;  [ 0.000 0.000  0.000 -0.000]
;  [ 0.000 0.000  0.000  0.000]]
; nil
; user=> (let [c (mmul m+ m)] (pm (sub c (transpose c))))
; [[ 0.000  0.000 0.000]
;  [-0.000  0.000 0.000]
;  [ 0.000 -0.000 0.000]]
