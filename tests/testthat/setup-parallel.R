## Test fixtures are deliberately small, and every dispatch in them lands under the
## `combat.min.cells` default. Left alone, the whole suite would run through the serial
## fallback and still report green: the backend matrix, the chunk-order tests and the
## equivalence tests would all be asserting on code that never forked.
##
## So the gate is off for the suite, and the two tests that exist to check the gate set
## their own value. Speed is not the point here, coverage is.
## Both thresholds, not just the first. `combat.min.disp.cells` defaults higher than
## `combat.min.cells`, so zeroing only one would leave the dispersion path serial
## throughout the suite: the same blind spot, one option later.
## combat.min.ls.cells is the limma least-squares gate and defaults far higher than the
## other two, because on lm.series' fast branch a split is four lm.fit calls where one would
## do. Left at its default the whole limma suite would run serially and still report green.
## EVERY gate, not the ones anybody remembered. Each new companion has brought its own
## threshold, and each time the suite kept passing while the tests it added ran through the
## serial fallback and asserted nothing about the split. The list below is the output of
##   grep -rho 'getOption("combat\.min[^)]*)' R/
## and adding a companion means adding its gate here in the same commit.
##   combat.min.cells        2e4   row splits generally
##   combat.min.disp.cells   3e4   the tagwise dispersion path
##   combat.min.ls.cells     6e6   limma's least-squares fast branch
##   combat.min.norm.cells   2e5   the edgeR normalisation column split
##   combat.min.order.cells  4e6   the rank and quantile stages inside normalisation
##   combat.min.dupcor.cells 5000  duplicateCorrelation's per-gene REML split
##   combat.min.glm.cells    1e5   glmFit's row split, which the shared 2e4 gate forked too eagerly
withr::local_options(combat.min.cells = 0, combat.min.disp.cells = 0,
                     combat.min.ls.cells = 0, combat.min.norm.cells = 0,
                     combat.min.order.cells = 0,
                     combat.min.dupcor.cells = 0, combat.min.glm.cells = 0,
                     .local_envir = testthat::teardown_env())
