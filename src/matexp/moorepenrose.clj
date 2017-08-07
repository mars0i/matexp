(ns matexp.moorepenrose
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.linear :as lin]
            [clojure.core.matrix.utils :as u]))


;; Questions:
;; How should tolerance be calculated? 
;; Complex matrices are not fully part of core.matrix at present, right?
;; So I should use transpose rather than hermitian-transpose?
;; Eventually this needs to use hermitian-transpose on complex matrices.
;; Could use doseq-indexed from Mikera's cljutils for the loop.
(defn pinv 
  "Moore-Penrose pseudoinverse of matrix m calculated using svd.  tolerance 
  defaults to 0 during calculation. Absolute values of singular values that 
  are less than or equal to tolerance will be treated as zero.  Note requires
  an implementation that includes svd."
  ([m] (pinv m 0.0)) ;; default tolerance should be larger--maybe calculated
  ([m tolerance]
   (let [treat-as-nonzero? (fn [x] (> (Math/abs x) tolerance))
         diag-pinv (fn [singular-vals rows cols] 
                     (let [smaller-dim (min rows cols)
                           sigma+ (m/ensure-mutable (m/zero-matrix rows cols))]
                       (u/doseq-indexed [x singular-vals i]
                          (when (treat-as-nonzero? x)
                            (m/mset! sigma+ i i (/ x))))
                       sigma+))
         [rows cols] (m/shape m)
         {:keys [U S V*]} (lin/svd m)
         _ (println (m/shape U) (m/shape m) (m/shape V*)) ; DEBUG
         S+ (diag-pinv S rows cols)]
     (println S+)
     (m/mmul (m/transpose V*) S+ (m/transpose U)))))
