# k-move-pricing datasets: the rule and list for R2 and R3

Status: sealed 2026-11-07, after its referees' reports and the second screens, for the maintainer's decision
(dec-B635). Companion to [k-move-pricing.md](k-move-pricing.md), whose step 2 guard and step 3 test read these
datasets (call 7, dec-B635, dec-B651).

The rule was sealed on 2026-11-07 by committing together, before any data file was loaded, this file, the per-entry
screen [k-move-pricing-screen.csv](k-move-pricing-screen.csv), the plan's sizing of R2 and R3 (ruled by method,
dec-B667, and read by a referee in two rounds), and, as copies in the folder k-move-pricing-seal beside this file,
verbatim but for the first sizing report's few non-ASCII symbols spelled in ASCII, listed with their checksums in
its MANIFEST.txt: the frame file, the snapshot checksums, the second screens (of the seven single-screened passers
and of stage 4's census entries, section 5) and the referees' reports. After the seal nothing here changes; an
error found later is reported in the study's report, not fixed by editing the rule. Nothing below rests on a data
file or an outcome column: only suite records, the TabArena curation record, repository metadata and published
papers were read.

## 1. What the datasets are for, and what that asks of them

Real outcomes enter pricing twice. R2 (12 datasets, sized in the plan's [Sizing the real-data
stages](k-move-pricing.md#sizing-the-real-data-stages)) carries step 2's real-data guard: each sampler's error at
equal running time against a long reference run of the same posterior ([The real-data
guard](k-move-pricing.md#the-real-data-guard)). R3 (every passer not in R2, the census read to its end first)
carries step 3's test of a trigger, if step 3 runs; R2's datasets join that test when the whole of step 3 is
committed before any R2 reading is opened (section 6).

Both readings compare two exact samplers of one posterior, dataset by dataset, against a reference of that same
posterior. No reading splits rows to estimate how a model generalizes, and no reading that frames a recommendation
scores an outcome the fit has not seen (the held-out score on the prediction rows is reported beside the guard,
framing nothing, and compares the two samplers on one fixed split). So a criterion matters here only if it bears on
a posterior comparison: whether the posterior is one a user meets, whether a choice of ours shaped it, whether two
datasets are one, and whether it can be fitted and used. The predictive-benchmark criteria the first rule copied
(independent rows, case-control sampling, leaks, duplicates, time-to-event status) protect a held-out score's claim
about new data, a claim this design no longer makes (section 3).

## 2. The frame

Four stages in a fixed, published order; stages 1 to 3 are screened in full, stage 4 only if they leave the pool
short (section 5).

| stage | source | version pinned | order within the stage |
|---|---|---|---|
| 1 | TabArena-v0.1, binary tasks (Erickson et al. 2025) | OpenML study 457, alias tabarena-v0.1, created 2025-05-20; 51 tasks, 30 binary | OpenML task id ascending |
| 2 | OpenML-CC18, binary tasks not seen and not in stage 1 (Bischl et al. 2021) | OpenML study 99, alias OpenML-CC18, created 2019-02-21; 72 tasks, 35 binary | OpenML task id ascending |
| 3 | binary entries of the TabArena curation record whose category is Tiny data, Image, Audio or Scientific Discovery and whose final decision is a "No" | tabarena.ai/dataset-curation as of 2026-11-07, 991 rows, sha256 ad6c9ec0afb440e4dca2bf90e27bb68e83ab59059c0290111e572ef9d239ded7 | record row order |
| 4 | seeded UCI census: every UCI entry listing Classification, at least 100 instances and Multivariate or Tabular, 387 entries | UCI API listing of 2026-11-07, 690 entries; permuted by `set.seed(20261107L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")` then `ids[sample.int(length(ids))]`, order file sha256 f8978e03450947106d2cacaf35fb5ed89aa2297889f5c511cc74a200f9b423d5 | the draw |

"Binary" in stages 1 and 2 is the suite's NumberOfClasses of 2. Stage 3's mechanics (one row per source under
several names; different files of one source entering apart and flagged) are as the frame file records them.
The frame file holds 538 rows (sha256 f6f3b0cb9780662fd9be71019da4baf37b536f84cbbbde89b521e242e4841055): the
151 binary entries of stages 1 to 3 in stage order (rows 1 to 151, seen datasets and earlier-stage copies among
them marked, 113 entering the screen), then the census's 387 entries in draw order (rows 152 to 538).

**Seen.** The 22 real prediction datasets of the k study (tictactoe, banknote, mushroom, spambase, ionosphere,
sonar, wdbc, climate, biopsy, magic, bank, adult, creditapproval, cleveland, pbc, kyphosis, pima, german,
transfusion, infert, birthwt, haberman), lindner, ACTG 175 and lalonde, and IHDP, ACIC 2016 and 2019 and Twins. A
suite copy of a seen dataset is seen (breast-w is biopsy, diabetes is pima, credit-g is german, and so on), and
so is a dataset from the same source, institution, collection or patients (Statlog heart and the Hungarian, VA
and Zurich heart-disease files beside cleveland; Census-Income beside adult; the corrected South German credit
beside german).

## 3. The screen

A doubt recorded by a suite's curators, an originating paper or either auditor fails the entry when it bears on a
kept criterion below or on the pre-fit check (section 4), counted under that criterion; a doubt on a dropped
matter (time order, grouped rows, leaks, duplicates, designed sampling, a task too easy, too small or outdated)
fails nothing, or the dropped criteria would return through it. The flow table counts each failure under its first
failing criterion in the order below.

| # | criterion | published source | why it bears on a posterior comparison |
|---|---|---|---|
| 1 | **Outcome.** The binary task as the originating study or donor defined it, with its predictors as the source lists them, their types included (the loader sets each predictor's type as the donor lists it, section 4, so a doubt on a suite file's types alone fails nothing); no merge, threshold or choice of target made by a later curator or by this study. A source with several candidate targets and no default fails. A suite's curated file is the file: its curators' row edits (duplicates removed, a corrupt row dropped) and value edits (missing values filled, training and test files combined) are allowed and recorded in the screen file, since they choose no target and no predictor. | CC18 excludes transformed targets (Bischl et al. 2021, section 6.1); AMLB excludes binarized tasks (Gijsbers et al. 2024); TabArena criterion 5 (Erickson et al. 2025) | a target or threshold of our choosing is a choice that could shape the posterior, even made blind |
| 2 | **Real data.** Not generated, simulated or augmented, with a named originating study or donor. | CC18 (a) and (d); TabArena criterion 4 | the maintainer's rule asks for interesting real problems; simulated data are what step 2's drawn truths already are |
| 3 | **Provenance.** Not from the same source, institution, collection or patients as a seen dataset; one dataset per source in the pool: of a source's entries, the earliest in frame order that passes every other criterion is kept (spect over spectf, which fails 5). | TabArena criterion 1; CC18 (b) | the held-out rule (dec-B530); two files of one source are one posterior twice |
| 4 | **Size and outcome share.** At least 100 rows and a listed minority share of at least 0.01. The share is listed by the suite for its file, by the donor's documentation, by the originating paper or by a repository record of the same file; a later user's count or a share derived from a summary statistic does not establish it, and a share not established fails. | CC18 (a), (c), (d) and TabArena criterion 6 set higher floors | keeps the small samples and rare outcomes the move is priced on; the floors are ours |
| 5 | **License and file.** The originating terms permit research use, read and recorded with their address. Terms are a license, or the originator's or originating repository's statement that the data are released for public or research use. "None", "N/A", a blank, a re-upload's license field alone, a citation request alone, or competition rules silent on reuse fail. A stable versioned file (a DOI or a versioned address), holding the described data, its sha256 recorded by the loader at its first download after the seal and checked at every later load: a file whose row count differs from the originating source's by no recorded edit fails (spectf, cylinder-bands). A source of several files enters as their concatenation (section 4). | TabArena criterion 8; CC18 "publicly licenced" | a precondition of using the data at all |
| 6 | **Feasibility.** At most 2,000 columns after expanding factors as bart() expands them (a two-level factor one column, a k-level factor k), from the listed feature metadata. | CC18 (b), fewer than 5,000 features | the reference runs' cost; the cap is ours |

Every kept criterion is one of population (which posteriors users meet, and whether a choice of ours shaped one)
or of feasibility. No criterion, kept or dropped, can favor one exact sampler over the other on a dataset: each is
a property of the data, and both samplers target the same posterior. What the criteria change is the population
of posteriors read and the share of references that pass.

No file holding any of these data, or derived from them, is committed or shipped; only summaries and readings
are (irish's terms allow research use only, and one rule serves every dataset).

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
  sampler except as near separation does, which is the effect measured. The draw by sealed score fixes R2's and
  R3's share of each rank block of the pool, not of each score range, so leaks that raise many datasets toward
  separation raise R2's near-separated share with the pool's, and with it the share of references that fail: a
  cost in power, not in validity. Scores of 0.995 or more stay in and are reported apart (section 6).
- *Duplicate rows on the covariates.* They change the design's effective size, not the posterior's validity.
- *Time-to-event status* (a death or recurrence flag at variable follow-up). The flag is a recorded binary
  outcome and its posterior is well defined; how well it predicts an event is not what is read.
- Never applied: the suites' "not too easy" filters, which select on the outcome and would remove the nearly
  separable data the move is priced on.

What the dropped criteria still cost: the pool is curated prediction tasks, not the observational and effect
data BART users often fit, and results hold for this frame only (Boulesteix, Wilson and Hapfelmeier 2017).

## 4. The covariate-only pre-fit check

After the seal and before any score, one loader, committed before it runs, reads each passer's file, records its
sha256 (checked at every later load), archives the file privately outside the tree, and removes the outcome column,
named by the suite's default target or the donor's documentation, immediately after parsing. It keeps the columns
the source lists as predictors: identifier columns, columns the donor marks as not predictors, and ignore
attributes are dropped (Mammographic Mass's BI-RADS assessment, which its donor marks non-predictive). It sets each
predictor's type, numeric or factor, as the donor lists it; where the donor lists none, or a type that is neither
(an ordinal, a date-time), it keeps the file's type, numbers numeric and text a factor. Applied to every passer
before the seal, this decides none (the screen file's types column): of the passers' donors only Mammographic
Mass's lists such a type, its density ordinal and coded as integers, so numeric; two suite files hold text columns
of a few levels (coil2000's age and contribution bands, seismic-bumps's hazard grade), factors either way; and of
the entries screened that failed, Occupancy Detection's date-time is text, a factor of about 20,560 levels, which
fails 6 after the 4 it fails first. Its committed report gives, from covariates only, the row and column counts
against the listing (an unexplained mismatch fails, as a doubt that the file is the described data), missing values
per column against the donor's, any value the suite filled, and columns after factor expansion.

- **The listing** is the source's stated row count, or the suite's for a suite file. Where a source states two,
  the screen file names the one taken: Parkinsons lists 197 instances in its metadata and 195 recordings in its
  description, and 195 is the listing. A suite file's rows are also set against the originating source's: a
  difference the curators' recorded edits give is explained (NATICUSdroid, deduplicated from 29,333 rows to
  7,491; qsar-biodeg, 1,054 against 1,055), and one they do not give fails criterion 5 (cylinder-bands, 540
  against the donor's 512).
- **Several files** are concatenated in the order the source lists them (spect: its training file of 80 rows,
  then its test file of 187).
- **Missing values** are kept, as a user's default fit keeps them: bart()'s `na.keepPredictors` routes a missing
  predictor down a learned side of each split, and rpart, which computes the separation score, uses surrogates.
  Only a column with no observed value is dropped. No row is dropped for a missing predictor, so the listing's
  counts and criterion 4's floors stand (OpenML lists a missing value in 75,244 of APSFailure's 76,000 rows, and
  OpenML's copy of hepatitis in 75 of its 155; deleting them would cut hepatitis below the floor). The hepatitis
  file used is PMLB's, whose missing values the check counts; a fill PMLB made is recorded, as a suite's value
  edit is. A row whose
  outcome is missing is dropped when the harness attaches the outcome, as `na.keepPredictors` drops it, before the
  rows are drawn; its count is in the harness's sealed record, and the minority share is not rechecked.

The pool is judged short on documentation before the seal and again after this check (section 5); the check, with
any stage-4 continuation, completes before any score (section 6).

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
| jm1 | "a PROMISE data set made publicly available in order to encourage repeatable, verifiable, refutable, and/or improvable predictive models", with the repository's acknowledgment guidelines | OpenML 1053 description, promise.site.uottawa.ca/SERepository | pass: the originating repository's statement of public release |
| analcatdata_boxing1, boxing2, creditscore | "may be used freely for scientific, educational and/or noncommercial purposes" (Simonoff 2003) | OpenML 448 description | pass on 5; fail 2, below |
| biomed | source unknown, no terms | OpenML 481 | fail |
| Pumpkin_Seeds | the creator's page gives a citation request and no terms; the Kaggle page it is listed from could not be read | muratkoklu.com/datasets | fail: a citation request alone; if terms are found and its share is established, it enters (rice, its source's later file, fails 4) |
| phoneme | a KEEL copyright notice, no terms; original from the ELENA project | sci2s.ugr.es/keel, dataset 105 | fail |
| electricity, SpeedDating, numerai28.6 | no terms on the named source pages | OpenML 151 and 40536; the Kaggle page states none | fail |
| mammography, oil_spill, sa-heart, prnn_crabs | terms not stated (doubts the first auditor recorded) | | fail |

Three settled points: analcatdata_lawsuit fails 2, the first auditor's unsettled doubt that its originating case
is unnamed, and boxing1, boxing2 and creditscore fail 2 on the same ground, since they come from the same book
with the same "Source: Unknown" record (one rule for one collection). PhishingWebsites fails 1 on the first
auditor's doubt that its outcome's encoding is not described, and molecular-biology_promoters fails 2 on the first
auditor's finding that its negatives are constructed. Parkinson_Multiple_Sound_Recording fails 1: its source lists
no target and offers the class and the UPDRS score.

**Flow**, each failure counted under its first failing criterion; stage 4's column is its first 80 draws, read
on one screen, its first four candidates on two:

| stage | 1 TabArena-v0.1 | 2 CC18 | 3 tier 2 | total, 1 to 3 | 4 census, 80 drawn |
|---|---|---|---|---|---|
| binary entries listed | 30 | 35 | 86 | 151 | 80 drawn, 43 not two-class |
| seen | 4 | 12 | 11 | 27 | 2 |
| already in an earlier stage | 0 | 4 | 7 | 11 | 11 |
| entering | 26 | 19 | 68 | 113 | 24 |
| fail 1, outcome | 0 | 3 | 22 | 25 | 6 |
| fail 2, real data | 9 | 3 | 12 | 24 | 6 |
| fail 3, provenance | 0 | 5 | 5 | 10 | 0 |
| fail 4, size or share | 0 | 0 | 11 | 11 | 3 |
| fail 5, license or file | 7 | 4 | 10 | 21 | 3 |
| fail 6, feasibility | 1 | 0 | 0 | 1 | 1 |
| pass | 9 | 4 | 8 | 21 | 5, 3 on two screens |
| of which at least 1,000 rows | 9 | 3 | 2 | 14 | 1, on one screen |

The per-entry table is [k-move-pricing-screen.csv](k-move-pricing-screen.csv): each of the 113 entering entries of
stages 1 to 3 and the 24 entering census entries read so far, in frame order, with its listing, its cell on each of
the six criteria (P or F, an unread or unsettled cell failing), its result as pass or its first failing criterion,
its second screen (done, or needed if step 3 runs), for a passer the predictor types read from its donor or
suite record (section 4), and the source of the deciding cell (for a passer, its license, its share's
source and any recorded edit, concatenation or missingness).

**The second screen** of the seven passers that rested on one screener read each on all six criteria without this
screen's cells. APSFailure, Internet-Advertisements and skin-segmentation pass. nomao fails 1, on the curators'
open doubt on its labels; clean1 fails 1, no source naming its task (read as Musk version 1, its donor's task is
per molecule); NHANES_age_prediction fails 1, its senior threshold being the subset preparer's, not the survey's;
cylinder-bands fails 5, its suite file holding 540 rows against the donor's 512 with no recorded edit to explain
them. Its readings of the rule, each now stated in sections 3 and 4 and applied to every entry (calls 11 to 15):
- a doubt on a dropped matter fails nothing. The first screen took kept-criterion cells only, and every doubt the
  curation record holds on the other passers is on a dropped matter (time order, grouped rows, a leak, small,
  outdated, an image domain), so nothing else changes;
- factors count as bart() expands them. The first screen counted a k-level factor as k - 1 columns, which only
  undercounts, and no passer comes near 2,000 columns, so nothing changes;
- a share is listed only by the suite, the donor, the originating paper or a repository record of the same file.
  Applied to every passer: the shares of spect, hepatitis and vote are now taken from their donors' files, of
  breast_cancer from UCI's record, of Parkinsons (0.246) and Diabetic_Retinopathy_Debrecen (0.469) from their
  originating papers, where the first screen had derived some from PMLB's imbalance statistic; rice fails 4, since
  no source of that kind states its classes (its originating paper was not reached and a later user's count does
  not count);
- a suite's value edits are recorded, not failed: Internet-Advertisements' filled continuous features and its
  one-row difference in class counts, APSFailure's combined files and its missing-value coding;
- predictors' types are the donor's, set by the loader (the screener read a departure as failing 1; call 15).
  cylinder-bands' type doubt is then resolved, and it fails on its file instead; no other entry changes.

**The 21 passers**, in frame order (rows as listed; the share is the listed minority share):

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
| 54 | 2 | ozone-level-8hr | 2,534 | 0.063 | environment |
| 62 | 2 | wilt | 4,839 | 0.054 | remote sensing |
| 64 | 2 | Internet-Advertisements | 3,279 | 0.140 | web |
| 70 | 3 | irish | 500 | 0.444 | education |
| 77 | 3 | skin-segmentation | 245,057 | 0.208 | imaging |
| 97 | 3 | spect | 267 | 0.206 | health |
| 101 | 3 | breast_cancer (Ljubljana) | 286 | 0.297 | health |
| 126 | 3 | hepatitis | 155 | 0.206 | health |
| 134 | 3 | Parkinsons | 195 | 0.246 | health |
| 136 | 3 | Diabetic_Retinopathy_Debrecen | 1,151 | 0.469 | health |
| 150 | 3 | vote | 435 | 0.386 | politics |

Fields: health 6; commerce, credit, finance and insurance 4; the sciences and engineering 8; others 3. Shares under
0.06: APSFailure, wilt.

**Against the orchestrator's count of 51 of 137, 27 with at least 1,000 rows.** That count passed every entry
whose license the second auditor had not read, included the census rows, did not apply the first auditor's
doubts or the reconciliation, and did not rerun the one-per-source rule. Counted that way over stages 1 to 3 the
figure is 43, 24 with at least 1,000 rows.

**The pool is short, so stage 4 runs.** 21 is under the bar of 24, and stage 4 runs in draw order until 24 pass or
the census ends, each census entry screened by section 3 by two auditors; no dataset is added by hand. It runs on
documentation before the seal, so the sealed list is whole, and its second screens are sealed with it. The census's
first 80 draws were read by the second auditor under the first rule and re-read under section 3, with the four
licenses that auditor left unread read from their UCI pages (all CC BY 4.0); a second screen then read the
candidates in draw order to the third pass:

| frame | dataset | rows | share | first screen | second screen |
|---|---|---|---|---|---|
| 152 | Toxicity (UCI 728) | 171 | 0.327 | pass | pass |
| 188 | Occupancy Detection (UCI 357) | 20,560 | | fails 4: share not established | fails 4, and 6 |
| 193 | Bone marrow transplant: children (UCI 565) | 187 | 0.450 | fails 4: share not established | pass, share settled |
| 200 | Mammographic Mass (UCI 161) | 961 | 0.463 | pass | pass, a class-count doubt recorded |
| 207 | Mesothelioma (UCI 351) | 324 | | fails 4: share not established | if step 3 runs: may settle the share |
| 211 | SECOM (UCI 179) | 1,567 | 0.066 | pass | if step 3 runs |
| 222 | Vertebral Column, two-class file (UCI 212) | 310 | 0.323 | pass | if step 3 runs |

The second screen's readings, set against the rule as written:
- Toxicity passes. Its page offers a second file of 13 features, the donors' feature selection on the outcome;
  the file is the one holding the 1,203 predictors the source lists.
- Occupancy Detection fails 4. No source lists the share of the file used: the donor's script comments give it per
  file only, and the screener recorded that doubt, which bears on a kept criterion. Its date predictor is listed as
  a date-time, neither numeric nor factor; under section 4 it keeps the file's text type, a factor of about 20,560
  levels, so it would fail 6 as well. Its outcome passes 1: an encoding is not a choice of target.
- Bone marrow transplant: children passes. Its target is the repository's default (survival_status), so
  criterion 1's no-default clause does not apply. The listing taken is 38 predictors: the 37 conditional
  attributes both donor papers give (Wrobel, Gudys and Sikora 2017, Table 1; Sikora, Wrobel and Gudys 2019) and
  survival_time, as the repository's roles enter it. survival_time is a leak for the status, a dropped matter that
  fails nothing. The other counts are recorded beside it, as Parkinsons's two counts are (section 4): the
  repository table's 36 features, survival_time among them; the donors' list of 39 entries, which adds
  survival_status; and the 2017 paper's text count of 42. A file holding other than 38 predictors fails the
  pre-fit check, and the census continues before the draw. Its share is the 2017 paper's, 55 percent of the 187
  censored, a share of one class, which criterion 4 accepts.
- Mammographic Mass passes. Its donor lists 516 benign and 445 malignant, its originating paper 515 and 446 (Elter,
  Schulz-Wendtland and Wittenberg 2007): one case, the rows agreeing at 961 and the share 0.463 or 0.464 by either.
  The doubt bears on no criterion's verdict and is recorded, not failed, as Internet-Advertisements's one-row class
  difference is. Its density is listed ordinal, kept as the file's integer codes under section 4.
- Two further readings apply to every entry and change none: the sha256 is the loader's to record (criterion 5),
  since no file is fetched before the seal; and a repository holding the creators' own deposit is the originating
  repository (criterion 5's terms), as every UCI license in section 5 was read.

Every other census entry among the first 80 draws fails on a fact or a recorded doubt that a second screen cannot
reverse (Record Linkage, for one, lists 20,931 matches among 5,749,132 pairs, a share of 0.004). The pool is 24,
14 with at least 1,000 rows, and the census stops there for R2.

**If step 3 runs**, the census is read on in draw order to its end (frame row 538) under both screens before any R3
file is loaded, and every further passer joins R3: R3 is every passer of the frame not in R2 (the plan's [Sizing
the real-data stages](k-move-pricing.md#sizing-the-real-data-stages)). These screens are committed before R3's
references run; they add nothing to R2.

**If the pre-fit check leaves the pool short again** (a file that fails its count or its checksum), stage 4
continues in draw order after the seal, under the sealed rule and two auditors, its screens committed before any
score; if the census ends short, the guard and the test run at the counts that pass.

## 6. Allocation, the draw and the order of reading

1. **First the pre-fit check.** After the seal, section 4's check runs on every passer, and stage 4 continues if
   the pool is then short, before any score is computed: the draw needs every passer's score.
2. **The separation score.** The harness computes each passer's score, the k study's 5-fold cross-validated rpart
   AUC (`separationScore` in [binary-k-mixing-need.R](../../benchmarks/R/binary-k-mixing-need.R): a
   classification tree, `xval = 0`, the AUC taken as at least 0.5), on its fitted rows (item 5), with the folds of
   the code below in place of that function's keyed seed, and keeps it sealed until step 3 starts.
3. **R2** is drawn from all passers. Sorted by score ascending, ties by frame order, passer of rank r is in block
   ceiling(n2 r / N); one is drawn from each block. n2 is 12, set by the plan's sizing
   ([Sizing the real-data stages](k-move-pricing.md#sizing-the-real-data-stages)); at 24 passers R2's blocks hold
   2 each. **R3** is every passer not in R2, with the census's further passers if step 3 runs (section 5), and is
   not drawn. No passer is held in reserve (item 6).
4. **The code**, in R 4.6, fixed now, with `frame` a dataset's frame order, `n` its rows after section 4, and `ord`
   the passers' names in score order:

   ```r
   rng <- function(seed) {
     set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
   }
   # rows, per dataset: m drawn, four fifths of them fitted, the rest prediction rows
   rng(20261109L + frame)
   m <- min(n, 5000L); rows <- sort(sample.int(n, m))
   fit <- sort(rows[sample.int(m, (4L * m) %/% 5L)]); pred <- setdiff(rows, fit)
   # the separation score's folds, on the fitted rows
   rng(20262109L + frame)
   folds <- sample(rep_len(1:5, length(fit)))
   # the allocation
   drawBlocks <- function(o, k) {
     b <- ceiling(k * seq_along(o) / length(o))
     o[vapply(seq_len(k), function(j) { i <- which(b == j); i[sample.int(length(i), 1L)] }, 1L)]
   }
   rng(20261108L)
   r2 <- drawBlocks(ord, 12L)
   r3 <- ord[!ord %in% r2]
   ```

   The allocation's seed and its block draw are the first rule's; the folds' seed is new.
5. **Rows.** For each dataset, the seeded simple random subsample of m = min(n, 5,000) rows above; four fifths are
   fitted (4,000 at most, the plan's cap) and the rest are prediction rows, whose outcomes stay sealed until the
   guard's reading is committed and are then read only for the held-out score, framing nothing.
6. **No reserves.** The pre-fit check (item 1) runs before the draw, and a shortfall it leaves is met by the
   census, so every drawn file has already loaded and matched its listing. A file that fails to reload after the
   draw is restored from the loader's private archive, or fetched again from its recorded address, and checked
   against its recorded sha256; one that cannot be is reported unresolved and not replaced. Nothing that fails or
   misbehaves after a score or a fit is replaced; a reference that fails its checks leaves the dataset unresolved,
   reported.
7. **Scores of 0.995 or more stay in**, in the top block, and every reading is also given without them: users fit
   such data, and the plan's rejection of them was a rule for drawing simulated truths.
8. **Order of reading.** The assignment is committed before step 2 runs. The whole of step 3 (the trigger itself,
   computed mechanically from step 2's simulated fits, the need-the-move rule, the refutation rule and
   R3's size) is committed before any reading of R2 is opened; R2's references and fits may run earlier,
   their summaries held closed. Step 2's report waits for that commit. Then R2's datasets, read for the guard, can
   also serve step 3's test, which reads every passer; any later amendment of step 3 drops R2 from it. The
   guard's reading is committed before the prediction rows' outcomes are unsealed.

## 7. What has no published model

Each is the study's own; the report says so.

| element | what it is | alternative |
|---|---|---|
| the dropped criteria | the argument of section 3 | the first rule's full predictive screen, which passed 7 of 231 |
| the 100-row and 0.01-share floors, the 2,000-column cap | keep small samples and rare outcomes; bound the reference's cost | the suites' floors of 500 rows and 0.05 share |
| the one-per-source tie rule | the earliest in frame order that passes every other criterion | the earliest in frame order outright, which fails spect behind spectf (20 pass) |
| the separation score and its blocks | the k study's rpart AUC, R2 and R3 stratified on it | stratify by size or outcome rate |
| the doubt-fails screen | a recorded doubt on a kept criterion or the pre-fit check fails | list doubted entries apart and read with them as a sensitivity analysis |
| the listed share | a share stated by the suite, the donor, the originating paper or a repository record of the same file | a share derived from a summary statistic or a later user's count |
| the donor's types | the loader sets each predictor's type as the donor lists it | fail a file whose types depart from the donor's listing |
| analcatdata's shared doubt | one doubt applied to the book's four files | screen each file alone (boxing1 and creditscore would pass) |
| the license standard | a statement of release for public or research use counts as terms; a citation request alone does not | an explicit license or terms of use only, or a citation request counted as terms |
| curated suite files | a suite's curated file is the file, its row edits recorded | fail any entry whose suite file differs from its source by row edits |
| no reserves | the pre-fit check before the draw, the census meeting its shortfall; each file archived at its first download, a file lost after the draw restored, re-fetched or reported unresolved | a fixed count of reserves, which the pre-fit check's unknown rate cannot size |
| a type neither numeric nor factor | the file's type kept, numbers numeric and text a factor | each such type read as a screener judges it |

## 8. Calls

Made here, each with its alternative. Counts are of stages 1 to 3 (21 pass, 14 with at least 1,000 rows; with the
census's three passers, 24 and 14); an alternative that removes a passer makes stage 4 read that much deeper
into the census.
1. The screen keeps outcome, real data, provenance, size, license and feasibility, and drops the rest (section 3).
   Alternative: keep the predictive screen, under which the pool held 7 passers, 2 with at least 1,000 rows,
   against the 24 the readings need.
2. R2 draws from passers of every size. The 1,000-row floor existed so a test fold held 200 rows, and the guard
   reads no test fold. Alternative: keep the floor (14 eligible, of which R2 at 12 would take 12); R2 would then
   reach larger real fits, beside step 2's simulated regimes, which are mostly under 1,000 rows.
3. The loader keeps the predictors as the source lists them (section 4). Alternative: drop only the outcome, which
   keeps columns a donor marks as not predictors, such as Mammographic Mass's BI-RADS assessment.
4. Licenses unread by the second auditor were read here (section 5), under one standard for every entry: a license,
   or the originator's or originating repository's statement that the data are released for public or research use,
   counts as terms, and a citation request alone does not; jm1 passes on PROMISE's statement and Pumpkin_Seeds
   fails. Alternatives: an explicit license or terms of use only, which removes jm1 (20 pass, 13 large); a citation
   request counted as terms, which admits Pumpkin_Seeds if its share is established (rice fails 4); or every
   unread license failed, which removes irish and jm1 (19 pass, 13 large).
5. The analcatdata files share lawsuit's doubt. Alternative: screen each alone.
6. The seeds are the ones the first rule fixed, unchanged: allocation 20261108, rows 20261109 plus frame order;
   the score's folds, new, 20262109 plus frame order (section 6).
7. Parkinson_Multiple_Sound_Recording fails 1, its source naming no target. Alternative: read the class as its
   target from the dataset's title and paper, which adds one passer of 1,040 rows.
8. Missing predictors are kept, as bart()'s default `na.keepPredictors` keeps them, and only an all-missing column
   is dropped (section 4): the posterior read is the one a user's default fit targets. Alternative: complete cases
   with the row floor re-applied, which leaves hepatitis about 80 rows and fails it (20 pass or fewer), cuts
   APSFailure to about 750 rows, and reads a posterior no default fit targets.
9. A suite's curated file, its row and value edits recorded, passes criterion 1 (NATICUSdroid's deduplication,
   qsar-biodeg's one row). Alternative: fail row edits, which removes those two (19 pass, 12 large) and any other
   curated file the pre-fit check finds.
10. Parkinsons's listing is 195, its description's count of recordings, not its metadata's 197. Alternative: 197,
    at which a file of 195 rows is an unexplained mismatch and Parkinsons fails at the pre-fit check (20 pass).
11. A recorded doubt fails only when it bears on a kept criterion or the pre-fit check (section 3). Alternative:
    every recorded doubt fails, which removes at least APSFailure, skin-segmentation, ozone-level-8hr,
    seismic-bumps, jm1, irish and ilpd (time order, grouped rows, a leak, outdated or a toy task) and the five
    small passers the curators called too small, leaving about 9, fewer if an image domain counts, and brings
    back through the curators' notes the criteria section 3 drops.
12. Factors count as bart() expands them (criterion 6). Alternative: full indicator coding, a two-level factor two
    columns, which fails Internet-Advertisements at 3,113 columns (20 pass, 13 large).
13. A share is established only by the suite, the donor, the originating paper or a repository record of the same
    file (criterion 4). Alternative: a share derived from a summary statistic or a later user's count, which
    passes rice on a later paper's 1,630 and 2,180 (22 pass, 15 large).
14. A suite's value edits are recorded, not failed (criterion 1). Alternative: an unrecorded value edit is a doubt
    that the file is the described data and fails at the pre-fit check, which removes Internet-Advertisements
    (20 pass, 13 large) and any PMLB file the check finds filled (hepatitis and breast_cancer list missing values).
15. The loader sets each predictor's type as the donor lists it (sections 3 and 4). Alternative, the second
    screener's: a file whose types depart from the donor's listing fails 1. That would need every passer's file
    types set against its donor's listing before the seal, which neither screen did, and fails an unknown number;
    cylinder-bands fails either way. The donor's types are the data as described, which is the posterior a user
    reading the documentation fits.
16. Stage 4 runs on documentation before the seal, since the pool is already short there, and its two screens are
    sealed with the list. Alternative: seal the 21 and run stage 4 after the seal under the sealed rule, which
    seals a list known to be short.
17. The census's second screen is read against the rule as written (section 5). Occupancy Detection fails 4, its
    share stated only per file; Bone marrow transplant: children passes on 38 predictors, the donors' 37
    conditional attributes and survival_time, its other counts recorded; Mammographic Mass passes, its one-case
    class-count difference recorded. Alternatives: Occupancy entered with its date-time as a number and the donor's
    per-file script comments read as a listed share (25 pass, 15 large); Bone marrow failed on its several counts
    as a recorded doubt (23 pass, the census read deeper); or Mammographic Mass failed under 5 on its class-count
    doubt, which by one rule for one case also fails Internet-Advertisements (22 pass, the census read deeper).
18. A predictor type the donor lists that is neither numeric nor factor keeps the file's type (section 4), applied
    to every passer and deciding none. Alternative: read each such type as its screener judged it, a date-time as a
    number and an ordinal as numeric, a choice of this study's made entry by entry, which the type rule exists to
    avoid.

For the maintainer, with dec-B635's decision, after the seal: the rule and the list as a whole, sealed once an
independent referee's second round on the sizing found nothing blocking. The sizes are not put to the maintainer:
as ruled (dec-B667), they are set by a sample-size justification in the plan's [Sizing the real-data
stages](k-move-pricing.md#sizing-the-real-data-stages), R2 at 12 for the guard's assurance of 0.8 and R3 every
other passer, the frame measured whole.

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
