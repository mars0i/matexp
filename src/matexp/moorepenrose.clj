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
   (let [treat-as-zero?  (fn [x] (< (Math/abs x) tolerance))
         diag-pinv (fn [s] (let [smaller-dim (apply min (mx/shape s))  ;; pinv for rectangular diagonal matrix
                                 newm (mx/transpose s)]
                             (dotimes [i smaller-dim] ; is it necessary to go through every i? if x is truly zero, are the rest as well?
                               (let [x (mx/mget newm i i)]
                                 (if (treat-as-zero? x)
                                   (when (not (zero? x)) (mx/mset! newm i i 0.0))
                                   (mx/mset! newm i i (/ x)))))))
         {:keys [U S V*]} (lin/svd m)
         S+ (diag-pinv S)]
     (mx/mmul (mx/transpose V*) S+ (mx/transpose U)))))
