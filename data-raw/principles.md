🌍 High‑level vision for a principled refactoring of fmrireg

Turn the pipeline into a collection of atomic, self‑describing objects that are operated on by small, single‑purpose functions.
Every object knows everything that is needed to use it; the orchestration layer never reverse‑engineers anything.

Concretely this means
	1.	Single source of truth
	•	Only one event‐sequence representation (event_seq → always matrix).
	•	Only one regressor type (regressor) – single_trial_regressor() becomes a thin alias.
	•	HRF objects always expose span, nbasis, scale_factor.
	2.	Encapsulation over reconstruction
	•	Contrast objects carry their own colind (design‑matrix columns) when they are built – no later searching.
	3.	Clear separation of concerns
	•	“Make data” (design matrices, regressors) is entirely separate from “fit model”.
	•	“How to split data” (runwise vs chunkwise) is delegated to pluggable fitters.
	4.	Functional composition
	•	Remove the monolithic branching inside fmri_lm_fit(); it becomes a 20‑line dispatcher that calls runwise_fit() or chunkwise_fit().

⸻

🔧 Medium‑level action plan

area	what to do	why / benefit
Event infrastructure	Unify on event_seq that always stores a matrix (even 1‑col). Keep event_term as a thin wrapper.	one code path for categorical / continuous / multi‑column events.
Cells & contrasts	cells() stays; contrast factories attachattr(weights,"colind") <- match(colnames(weights), colnames(design_matrix(event_model)))	eliminates fragile lookup in fmri_lm_fit(); contrasts are now self‑contained.
HRF/Convolution	Delete convolve_block() and expose it through gen_hrf_blocked() only; pre‑compute Toeplitz/FFT helpers inside HRF object if speed is needed.	one obvious entry‑point; no duplicate code.
Regressors	regressor() takes care of single‑trial vs multi‑trial, amplitudes, durations, summation.single_trial_regressor() ← alias.	simplifies evaluate methods and plotting.
Model fitting	Split current fmri_lm_fit() intorunwise_fit() and chunkwise_fit(). Dispatcher does:fitter <- switch(strategy, ...) ; result <- fitter(...)	removes 80 + lines of nested conditionals; easier to test & extend.
Fast vs robust paths	Fitters decide: if all contrasts have colind → call C++ fast routine; else fallback to R implementation.	no global use_fast_path flag, fewer incompatible option combos.
Error handling / assertions	Keep assertions close to object construction; downstream code trusts its inputs.	defensive programming once, not everywhere.
Naming / clarity	functions describe action (evaluate_regressor_fast, runwise_fit) ; objects describe data (hrf, contrast, regressor).	discoverability & code navigation.


⸻

✨ Guiding principles
	•	Object first, helper later – build rich, immutable objects; helpers only use them.
	•	One way to do one thing – no duplicated convolution or event handling code.
	•	Fail fast & locally – validate in constructors, not in deep call stacks.
	•	Small surface, pluggable back‑ends – new fitters or HRF engines can be added without touching high‑level API.
	•	Readable > clever – explicit helpers beat magic branching; favour clarity even if you lose a micro‑second.

Implementing the first two bullet points (contrast encapsulation & fitter split) already eliminates the current crash and slices fmri_lm_fit() to a fraction of its size; the remaining steps round out the refactor into a clean, maintainable core.