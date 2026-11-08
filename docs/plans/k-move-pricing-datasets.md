# k-move-pricing datasets: the rule and list for R2 and R3

Status: draft for the referee, not sealed (2026-11-07). Companion to [k-move-pricing.md](k-move-pricing.md), whose
step 2 guard and step 3 test read these datasets (call 7, dec-B635, dec-B651).

The rule is sealed when this file, the frame file and the screen file it names, and the snapshot checksums are
committed together with the referee's report, before any data file is loaded. After the seal nothing here
changes; an error found later is reported in the study's report, not fixed by editing the rule. Nothing below
rests on a data file or an outcome column: only suite records, the TabArena curation record, repository
metadata and published papers were read.

## 1. What the datasets are for, and what that asks of them

Real outcomes enter pricing twice. R2 (6 datasets) carries step 2's real-data guard: each sampler's error at
equal running time against a long reference run of the same posterior ([The real-data
guard](k-move-pricing.md#the-real-data-guard)). R3 (12 datasets) carries step 3's test of a trigger, if step 3
runs; R2's datasets join that test when the trigger is sealed before R2 is read (section 6).

Both readings compare two exact samplers of one posterior, dataset by dataset, against a reference of that same
posterior. No reading splits rows to estimate how a model generalizes, and no reading scores an outcome the fit
has not seen. So a criterion matters here only if it bears on a posterior comparison: whether the posterior is
one a user meets, whether a choice of ours shaped it, whether two datasets are one, and whether it can be fitted
and used. The predictive-benchmark criteria the first rule copied (independent rows, case-control sampling,
leaks, duplicates, time-to-event status) protect a held-out score's claim about new data, a claim this design no
longer makes (section 3).

## 2. The frame

Four stages in a fixed, published order; stages 1 to 3 are screened in full, stage 4 only if they leave the pool
short (section 5).

| stage | source | version pinned | order within the stage |
|---|---|---|---|
| 1 | TabArena-v0.1, binary tasks (Erickson et al. 2025) | OpenML study 457, alias tabarena-v0.1, created 2025-05-20; 51 tasks, 30 binary | OpenML task id ascending |
| 2 | OpenML-CC18, binary tasks not seen and not in stage 1 (Bischl et al. 2021) | OpenML study 99, alias OpenML-CC18, created 2019-02-21; 72 tasks, 35 binary | OpenML task id ascending |
| 3 | binary entries of the TabArena curation record whose category is Tiny data, Image, Audio or Scientific Discovery and whose final decision is a "No" | tabarena.ai/dataset-curation as of 2026-11-07, 991 rows, sha256 ad6c9ec0afb440e4dca2bf90e27bb68e83ab59059c0290111e572ef9d239ded7 | record row order |
| 4 | seeded UCI census: every UCI entry listing Classification, at least 100 instances and Multivariate or Tabular, 387 entries | UCI API listing of 2026-11-07, 690 entries; permuted by `set.seed(20261107L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")` then `ids[sample.int(length(ids))]`, order file sha256 f8978e03...b423d5 | the draw |

"Binary" in stages 1 and 2 is the suite's NumberOfClasses of 2. Stage 3's mechanics (one row per source under
several names; different files of one source entering apart and flagged) are as the frame file records them.
The frame file holds 538 rows (sha256 f6f3b0cb9780662fd9be71019da4baf37b536f84cbbbde89b521e242e4841055):
113 entries from stages 1 to 3 after seen datasets and earlier-stage copies are removed, then the census in draw
order.

**Seen.** The 22 real prediction datasets of the k study (tictactoe, banknote, mushroom, spambase, ionosphere,
sonar, wdbc, climate, biopsy, magic, bank, adult, creditapproval, cleveland, pbc, kyphosis, pima, german,
transfusion, infert, birthwt, haberman), lindner, ACTG 175 and lalonde, and IHDP, ACIC 2016 and 2019 and Twins. A
suite copy of a seen dataset is seen (breast-w is biopsy, diabetes is pima, credit-g is german, and so on), and
so is a dataset from the same source, institution, collection or patients (Statlog heart and the Hungarian, VA
and Zurich heart-disease files beside cleveland; Census-Income beside adult; the corrected South German credit
beside german).

## 3. The screen

A doubt recorded by a suite's curators, an originating paper or either auditor fails the entry; the flow table
counts each failure under its first failing criterion in the order below.

| # | criterion | published source | why it bears on a posterior comparison |
|---|---|---|---|
| 1 | **Outcome.** The binary task as the originating study or donor defined it, with its predictors as the source lists them; no merge, threshold or choice of target made by a later curator or by this study. A source with several candidate targets and no default fails. | CC18 excludes transformed targets (Bischl et al. 2021, section 6.1); AMLB excludes binarized tasks (Gijsbers et al. 2024); TabArena criterion 5 (Erickson et al. 2025) | a target or threshold of our choosing is a choice that could shape the posterior, even made blind |
| 2 | **Real data.** Not generated, simulated or augmented, with a named originating study or donor. | CC18 (a) and (d); TabArena criterion 4 | the maintainer's rule asks for interesting real problems; simulated data are what step 2's drawn truths already are |
| 3 | **Provenance.** Not from the same source, institution, collection or patients as a seen dataset; one dataset per source in the pool, the earliest in frame order kept. | TabArena criterion 1; CC18 (b) | the held-out rule (dec-B530); two files of one source are one posterior twice |
| 4 | **Size and outcome share.** At least 100 rows and a listed minority share of at least 0.01. | CC18 (a), (c), (d) and TabArena criterion 6 set higher floors | keeps the small samples and rare outcomes the move is priced on; the floors are ours |
| 5 | **License and file.** The originating terms permit research use, read and recorded with their address; "None", "N/A", a blank, a re-upload's label or competition rules silent on reuse fail. A stable versioned file with a recorded sha256. | TabArena criterion 8; CC18 "publicly licenced" | a precondition of using the data at all |
| 6 | **Feasibility.** At most 2,000 columns after expanding factors, from the listed feature metadata. | CC18 (b), fewer than 5,000 features | the reference runs' cost; the cap is ours |

**Dropped, each with why it cannot bias the paired comparison.** Every reading compares two samplers whose
stationary distribution is the same posterior, on the rows as fitted, against a reference of that posterior.
- *Independent rows* (no time order, no repeated subjects, no rows grouped by entity). Dependence makes the model
  misspecified, but the posterior under it is still the one a user who fits these rows gets, and the one both
  samplers and the reference target. Nothing is split to estimate error on new data, which is what dependence
  would bias. Rows that repeat an entity concentrate the posterior as a larger sample would; the stratified draw
  by separation score (section 6) places such a dataset by where its posterior sits, not by why.
- *Case-control sampling and designed class balance.* They shift the intercept's meaning and the outcome rate,
  not the validity of the posterior given the sample; users fit logistic models to such samples. The outcome rate
  is a property the score blocks spread, not a defect.
- *Leaks* (post-outcome columns, preprocessing fitted on all rows). A leaking column makes the data nearly
  separable, and its posterior is one of the near-separated kind the move is priced on. It cannot favor either
  sampler except as near separation does, which is the effect measured, and the draw by sealed score fixes how
  many datasets come from each part of the score range, so leaks cannot raise R2's or R3's share of
  near-separated data. Scores of 0.995 or more stay in and are reported apart (section 6).
- *Duplicate rows on the covariates.* They change the design's effective size, not the posterior's validity.
- *Time-to-event status* (a death or recurrence flag at variable follow-up). The flag is a recorded binary
  outcome and its posterior is well defined; how well it predicts an event is not what is read.
- Never applied: the suites' "not too easy" filters, which select on the outcome and would remove the nearly
  separable data the move is priced on.

What the dropped criteria still cost: the pool is curated prediction tasks, not the observational and effect
data BART users often fit, and results hold for this frame only (Boulesteix, Wilson and Hapfelmeier 2017).

## 4. The covariate-only pre-fit check

After the seal and before any score, one loader, committed before it runs, reads each passer's file, verifies its
sha256, and removes the outcome column, named by the suite's default target or the donor's documentation,
immediately after parsing. It keeps the columns the source lists as predictors: identifier columns, columns the
donor marks as not predictors, and ignore attributes are dropped (NHANES's respondent age, the source of its
age-group target, is listed by the donor as "Other"). Its committed report gives, from covariates only, the row
and column counts against the listing (an unexplained mismatch fails, as a doubt that the file is the described
data), missing values per column and columns after factor expansion. Missing values are handled on covariates
only: a column more than half missing is dropped, then rows with any missing covariate. The pool is judged short
or not after this check (section 5).

## 5. The re-screen and the list

The 137 rows screened from documentation by the second auditor (113 from stages 1 to 3, 24 of the census's first
80 draws) were re-read under section 3. Cells already settled under the first rule carry over: the reconciliation
of the two audits and the referee's corrections stand, and a doubt either auditor recorded on a kept criterion
still fails. Licenses the second auditor left unread, because the entry failed elsewhere, were read here from the
originating pages:

| entry | terms read | address | result |
|---|---|---|---|
| irish | "may be used solely as an example for research on the development of statistical methods" | OpenML 451 description, owners Greaney and Kelleghan | pass |
| clean1 (musk, version 1) | CC BY 4.0 | archive.ics.uci.edu/dataset/74 | pass |
| jm1 | the PROMISE acknowledgment guidelines for published work | OpenML 1053 description, promise.site.uottawa.ca/SERepository | pass, flagged for the referee |
| analcatdata_boxing1, boxing2, creditscore | "may be used freely for scientific, educational and/or noncommercial purposes" (Simonoff 2003) | OpenML 448 description | pass on 5; fail 2, below |
| biomed | source unknown, no terms | OpenML 481 | fail |
| Pumpkin_Seeds | the creator's page asks for a citation and states no terms | muratkoklu.com/datasets | fail; if terms are found it enters and rice, its source's later file, leaves |
| phoneme | a KEEL copyright notice, no terms; original from the ELENA project | sci2s.ugr.es/keel, dataset 105 | fail |
| electricity, SpeedDating, numerai28.6 | no terms on the named source pages | OpenML 151 and 40536; the Kaggle page states none | fail |
| mammography, oil_spill, sa-heart, prnn_crabs | terms not stated (doubts the first auditor recorded) | | fail |

Three settled points: analcatdata_lawsuit fails 2, the first auditor's unsettled doubt that its originating case
is unnamed, and boxing1, boxing2 and creditscore fail 2 on the same ground, since they come from the same book
with the same "Source: Unknown" record (one rule for one collection). PhishingWebsites fails 1 on the first
auditor's doubt that its outcome's encoding is not described, and molecular-biology_promoters fails 2 on the first
auditor's finding that its negatives are constructed. Parkinson_Multiple_Sound_Recording fails 1: its source lists
no target and offers the class and the UPDRS score.

**Flow**, each failure counted under its first failing criterion:

| stage | 1 TabArena-v0.1 | 2 CC18 | 3 tier 2 | total |
|---|---|---|---|---|
| binary entries listed | 30 | 35 | 86 | 151 |
| seen | 4 | 12 | 11 | 27 |
| already in an earlier stage | 0 | 4 | 7 | 11 |
| entering | 26 | 19 | 68 | 113 |
| fail 1, outcome | 0 | 2 | 20 | 22 |
| fail 2, real data | 9 | 3 | 12 | 24 |
| fail 3, provenance | 0 | 5 | 5 | 10 |
| fail 4, size or share | 0 | 0 | 10 | 10 |
| fail 5, license or file | 7 | 3 | 10 | 20 |
| fail 6, feasibility | 1 | 0 | 0 | 1 |
| pass | 9 | 6 | 11 | 26 |
| of which at least 1,000 rows | 9 | 4 | 4 | 17 |

The pool is not short (24 passers is the bar), so stage 4 does not run, and the census rows read under the first
rule (Toxicity, vertebral column, mammographic mass and the others) are not in the pool. The per-entry table of
the re-screen goes into the sealed screen file.

**The 26 passers**, in frame order (rows as listed; the share is the listed minority share, blank where no source
read states it):

| frame | stage | dataset | rows | share | field |
|---|---|---|---|---|---|
| 2 | 1 | APSFailure | 76,000 | 0.018 | vehicle maintenance |
| 8 | 1 | coil2000_insurance_policies | 9,822 | 0.060 | insurance |
| 10 | 1 | credit_card_clients_default | 30,000 | 0.221 | credit |
| 24 | 1 | NATICUSdroid | 7,491 | 0.350 | software security |
| 25 | 1 | online_shoppers_intention | 12,330 | 0.155 | commerce |
| 26 | 1 | polish_companies_bankruptcy | 5,910 | 0.069 | finance |
| 27 | 1 | qsar-biodeg | 1,054 | 0.337 | chemistry |
| 28 | 1 | seismic-bumps | 2,584 | 0.066 | mining |
| 30 | 1 | jm1 | 10,885 | 0.193 | software |
| 51 | 2 | ilpd | 583 | 0.286 | health |
| 53 | 2 | nomao | 34,465 | 0.286 | record linkage |
| 54 | 2 | ozone-level-8hr | 2,534 | 0.063 | environment |
| 58 | 2 | cylinder-bands | 540 | 0.422 | manufacturing |
| 62 | 2 | wilt | 4,839 | 0.054 | remote sensing |
| 64 | 2 | Internet-Advertisements | 3,279 | 0.140 | web |
| 70 | 3 | irish | 500 | 0.444 | education |
| 77 | 3 | skin-segmentation | 245,057 | 0.208 | imaging |
| 97 | 3 | spect | 267 | 0.206 | health |
| 101 | 3 | breast_cancer (Ljubljana) | 286 | 0.297 | health |
| 118 | 3 | clean1 (musk, version 1) | 476 | | chemistry |
| 126 | 3 | hepatitis | 155 | 0.206 | health |
| 134 | 3 | Parkinsons | 197 | | health |
| 136 | 3 | Diabetic_Retinopathy_Debrecen | 1,151 | 0.469 | health |
| 138 | 3 | NHANES_age_prediction | 6,287 | | health survey |
| 143 | 3 | rice_cammeo_and_osmancik | 3,810 | | agriculture |
| 150 | 3 | vote | 435 | 0.386 | politics |

Fields: health 7; commerce, credit, finance and insurance 4; the sciences and engineering 11; others 4.
Shares under 0.06: APSFailure, wilt. Seven passers (APSFailure, nomao, cylinder-bands, Internet-Advertisements,
skin-segmentation, clean1, NHANES) rest on the second auditor's screen and this re-reading of their source pages;
the first auditor did not reach them, and the referee is asked to check them as the second screen.

**Against the orchestrator's count of 51 of 137, 27 with at least 1,000 rows.** That count passed every entry
whose license the second auditor had not read, included the census rows, did not apply the first auditor's
doubts or the reconciliation, and did not rerun the one-per-source rule. Counted that way over stages 1 to 3 the
figure is 43, 24 with at least 1,000 rows.

**If the pool is short.** If the pre-fit check leaves fewer than 24 passers, stage 4 runs in draw order, each
census entry screened by section 3 by two auditors and then checked, until 24 pass or the census ends; no dataset
is added by hand, and if the pool is still short the guard and the test run at the counts that pass.

## 6. Allocation, the draw and the order of reading

1. The harness computes each passer's separation score (a 5-fold cross-validated rpart AUC, the k study's) on its
   fitted rows after section 4, and keeps it sealed until step 3 starts.
2. **R2** is drawn from all passers. Sorted by score ascending, ties by frame order, passer of rank r is in block
   ceiling(6 r / N); one is drawn from each block.
3. **R3** is drawn from the passers not in R2: block ceiling(12 r / M) by the same order, one from each.
4. **The draw**, in R 4.6, fixed now, with `ord` the passers in score order:

   ```r
   set.seed(20261108L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
   N <- length(ord); b2 <- ceiling(6 * seq_len(N) / N)
   r2 <- ord[vapply(1:6, function(b) { i <- which(b2 == b); i[sample.int(length(i), 1L)] }, 1L)]
   rest <- ord[!ord %in% r2]; M <- length(rest); b3 <- ceiling(12 * seq_len(M) / M)
   r3 <- rest[vapply(1:12, function(b) { i <- which(b3 == b); i[sample.int(length(i), 1L)] }, 1L)]
   ```

   At 26 passers R2's blocks hold 4 or 5 and R3's 1 or 2 of the remaining 20; 8 are reserves.
5. **Rows.** For each dataset, one seeded simple random subsample of m = min(n, 5,000) rows, seed 20261109 plus its
   frame order; a seeded 80 percent are fitted (at most 4,000, the plan's cap) and the other 20 percent are
   prediction rows, whose outcomes are dropped at load and never read.
6. **Reserves** replace only a dataset whose file fails the pre-fit check when the harness reloads it, drawn from
   the same block. Nothing that fails or misbehaves after a score or a fit is replaced; a reference that fails its
   checks leaves the dataset unresolved, reported.
7. **Scores of 0.995 or more stay in**, in the top block, and every reading is also given without them: users fit
   such data, and the plan's rejection of them was a rule for drawing simulated truths.
8. **Order of reading.** The assignment is committed before step 2 runs. Step 3's trigger form is computed from
   step 2's simulated fits and committed before any reading of R2 is opened; R2's references and fits may run
   earlier, their summaries held closed. Then R2's datasets, read for the guard, can also serve step 3's test,
   which reads 18 datasets in place of 12.

## 7. What has no published model

Each is the study's own; the report says so.

| element | what it is | alternative |
|---|---|---|
| the dropped criteria | the argument of section 3 | the first rule's full predictive screen, which passed 7 of 231 |
| the 100-row and 0.01-share floors, the 2,000-column cap | keep small samples and rare outcomes; bound the reference's cost | the suites' floors of 500 rows and 0.05 share |
| the one-per-source tie rule | earliest in frame order | a seeded draw among a source's files |
| the separation score and its blocks | the k study's rpart AUC, R2 and R3 stratified on it | stratify by size or outcome rate |
| the doubt-fails screen | any recorded doubt fails | list doubted entries apart and read with them as a sensitivity analysis |
| analcatdata's shared doubt | one doubt applied to the book's four files | screen each file alone (boxing1 and creditscore would pass) |

## 8. Calls

Made here, each with its alternative:
1. The screen keeps outcome, real data, provenance, size, license and feasibility, and drops the rest (section 3).
   Alternative: keep the predictive screen, under which the pool held 7 passers, 2 with at least 1,000 rows,
   against the 24 the readings need.
2. R2 draws from passers of every size. The 1,000-row floor existed so a test fold held 200 rows, and the guard
   reads no test fold. Alternative: keep the floor (17 eligible); R2 would then reach larger real fits, beside step
   2's simulated regimes, which are mostly under 1,000 rows.
3. The loader keeps the predictors as the source lists them (section 4). Alternative: drop only the outcome, which
   leaves NHANES's target computable from a covariate.
4. Licenses unread by the second auditor were read here (section 5), jm1 passed on the PROMISE acknowledgment
   guidelines. Alternative: fail every unread license, which removes irish, clean1 and jm1 (23 pass, 16 large) and
   makes the pool short.
5. The analcatdata files share lawsuit's doubt. Alternative: screen each alone.
6. The seeds are the ones the first rule fixed, unchanged: allocation 20261108, rows 20261109 plus frame order.
7. Parkinson_Multiple_Sound_Recording fails 1, its source naming no target. Alternative: read the class as its
   target from the dataset's title and paper, which adds one passer of 1,040 rows.

For the maintainer, with dec-B635's decision: the rule and the list as a whole, after the referee; and the first
rule's call on R3's size (12, or every passer not in R2), which the shared reading of section 6 makes 18 or 26.

## References

Read for this rule:
- Bischl, B., Casalicchio, G., Feurer, M., et al. (2021). OpenML benchmarking suites. NeurIPS 2021 Datasets and
  Benchmarks; arXiv:1708.03731v3.
- Boulesteix, A.-L., Wilson, R. and Hapfelmeier, A. (2017). Towards evidence-based computational statistics.
  BMC Medical Research Methodology 17, 138.
- Erickson, N., Purucker, L., Tschalzev, A., et al. (2025). TabArena: a living benchmark for machine learning on
  tabular data. NeurIPS 2025 Datasets and Benchmarks; arXiv:2506.16791v4; and its curation record.
- Gijsbers, P., Bueno, M. L. P., Coors, S., et al. (2024). AMLB: an AutoML benchmark. JMLR 25.
- Magnusson, M., Torgander, J., Burkner, P.-C., Zhang, L., Carpenter, B. and Vehtari, A. (2024). posteriordb:
  testing, benchmarking and developing Bayesian inference algorithms. arXiv:2407.04967v1.
- The repository pages named in section 5.
