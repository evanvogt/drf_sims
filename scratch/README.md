# scratch

Unmaintained exploratory code. Kept because it records how something was worked
out, not because it runs.

| file | was for |
|---|---|
| `scratch_CSH.R` | cause-specific hazard data generation, before `competing_risk/surv_dgm.R` settled |
| `scratch_SL.R` | SuperLearner behaviour on small samples |
| `scratch_combine_results.R` | an earlier results-collection approach, superseded by `R/pipeline.R` |

**Nothing here is sourced by any study.** These files define functions that
shadow real ones (`hazard()` also exists in `competing_risk/surv_dgm.R` with a
different signature), so do not source them alongside study code.

They are excluded from the de-duplication work deliberately — consolidating
scratch code would only make it look authoritative.

`.gitignore` ignores `scratch/*test*`, which is where genuinely throwaway test
files should go. That rule used to be repo-wide `*test*.*` and was silently
excluding `crossfitting/cf_testing.R`, the verification harness the crossfitting
README tells you to run before submitting anything.
