(ns matexp.moorepenrose
  (:require [clojure.core.matrix :as mx]
            [clojure.core.matrix.linear :as lin]))


;; Questions:
;; How should tolerance be calculated? 
;; Complex matrices are not fully part of core.matrix at present, right?
;; So I should use transpose rather than hermitian-transpose?
;; Eventually this needs to use hermitian-transpose on complex matrices.
(defn pinv 
  "Moore-Penrose pseudoinverse of matrix m, calculated using svd.  tolerance 
  defaults to [FIXME]; during calculation, values below tolerance in the sigma
  matrix produced by svd are treated as zero."
  ([m] (pinv m 0.0)) ;; default tolerance should be larger--maybe calculated
  ([m tolerance]
   (let [treat-as-nonzero?  (fn [x] (> (Math/abs x) tolerance))
         diag-pinv (fn [singular-vals sigma-rows sigma-cols] 
                     (let [smaller-dim (min rows cols)
                           sigma+ (ensure-mutable (zero-matrix cols rows))]
                       (doseq [i (range smaller-dim)
                               :let [x (mx/mget singular-vals i)]
                               :while (treat-as-nonzero? x)]
                         (mx/mset! sigma+ i i (/ x)))))
         {:keys [U S V*]} (lin/svd m)
         rows (second (shape U)) ; FIXME
         cols (first (shape V*)) ; FIXME
         S+ (diag-pinv S rows cols)]
     (mx/mmul (mx/transpose V*) S+ (mx/transpose U)))))
