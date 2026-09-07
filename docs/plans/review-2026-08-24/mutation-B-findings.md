# Mutation testing, leg B: the C++ component suite and the engine

Snapshot at b102e17c, 2026-08-24; counts below are as of that commit.

At HEAD (8e8a63ad): tests/cpp/test_facade.cpp now drives every SamplerBase
virtual through the base class (about 1,200 lines) and covers the gaps
BLOCKER #2 names - numSavedDraws vs capacity, setResponse's updateScale,
savedSlotForDraw, savedTree, and setForestWeights all have dedicated checks
in that file. The "Ladder audit" section's premise - `default:` arms on the
ResponseFamily switches - no longer holds: no `default:` label remains in
src/bartcore/*.hpp or src/R_interface_bartcore.cpp. Other findings below
have not been re-audited since this snapshot.

Read-only at b102e17c.  Work in a `git archive HEAD` copy under
a private scratch prefix with the checkout's src/*.a and generated
config headers copied in; the checkout's src/ was never touched, and a final
`diff -r` confirmed every staged header restored byte for byte.  Long-form
evidence: mutation-B-evidence.md, same directory.  SCOPE the ten tests/cpp
files touched since 03b97db7 (1309 of the suite's ~1520 check sites) and the
nine engine headers they exercise.  METHOD 80 SEMANTIC mutations, one header
each: apply, `make -j8`, run the WHOLE suite plain (~16s clean, "all tests
passed"), record the first check() that fired, restore and touch; mean cycle
52s.  Three extra legs, each once: a 7th-enumerator compile probe, an
ASAN/UBSAN pristine run, and a -DNDEBUG re-run of the two mutations whose catch
could have needed a live assert().  RESULT 80 planted / 63 caught (47 by a
named check, 16 only by a crash) / 17 missed, 2 of those equivalent mutants -
15 genuine gaps.  Per header, planted(missed): scan 10(0) tree 7(0) grow 9(1)
combiner 8(1) chain 10(2) sampler 9(2) moves 12(3) model 8(3) facade 7(5).

## Mutation table
id   header:line       mutation                                                   verdict / caught by
S1   [src/bartcore/scan.hpp:175](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L175)      swap the two missing-routed entries (direction convention  caught every entry scores its cut with the missing rows in the ch
S2   [src/bartcore/scan.hpp:159](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L159)      drop the paired occupancy sentinel on entry[1] (adjudicati caught crash SIGSEGV, no check fired
S3   [src/bartcore/scan.hpp:146](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L146)      off-by-one in the prefix scan: codes 0..cut+1 go left      caught crash SIGSEGV, no check fired
S4   [src/bartcore/scan.hpp:158](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L158)      occupancy compare relaxed to < 0 (zero-weight side becomes caught a zero-weight-only side is vetoed, member count notwithsta
S5   [src/bartcore/scan.hpp:243](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L243)      exact partition decode off-by-one (candidate bit shift)    caught no induced partition is emitted twice
S6   [src/bartcore/scan.hpp:129](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L129)      missing rows routed into ordinal bin 0 instead of the miss caught a node holding missing rows emits both directions of every
S7   [src/bartcore/scan.hpp:178](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L178)      emitted count always numCuts (doubled layout mis-read as p caught a node holding missing rows emits both directions of every
S8   [src/bartcore/scan.hpp:304](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L304)      categorical histogram: stale `seen` cache never cleared    caught twelve levels run past the cap
G1   [src/bartcore/grow.hpp:292](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L292)      drop the -log 2 the doubled layout owes (adjudication P2)  caught realized root rules match the exact law over the enumerate
G2   [src/bartcore/grow.hpp:299](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L299)      candidate cut index not de-interleaved under the doubled l caught every deeply grown ordinal rule sits inside its ancestor i
G3   [src/bartcore/grow.hpp:302](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L302)      candidate's named missing direction inverted               caught realized root rules match the exact law over the enumerate
G4   [src/bartcore/grow.hpp:120](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L120)      inline absent-category coin dropped (absent categories pin caught an absent reachable category is drawn, not pinned to one s
G5   [src/bartcore/grow.hpp:222](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L222)      split-variable normalizer uses numPredictors, not numAvail MISSED -
G6   [src/bartcore/grow.hpp:341](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L341)      installed rule takes the OPPOSITE missing direction from t caught realized root rules match the exact law over the enumerate
G7   [src/bartcore/grow.hpp:284](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L284)      cut buffer sized numCuts, not 2*numCuts (layout/allocation caught crash SIGBUS, no check fired
G8   [src/bartcore/grow.hpp:201](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L201)      no-split candidate loses its (1 - growth) prior factor     caught realized root rules match the exact law over the enumerate
M1   [src/bartcore/moves.hpp:120](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L120)     veto rank compare >= on the current side (every proposal w caught dart concentrates splits on signal variables
M2   [src/bartcore/moves.hpp:128](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L128)     drop the both-infeasible NaN guard in resolveVetoRank      caught two sentinels at equal rank reject rather than differencin
M3   [src/bartcore/moves.hpp:281](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L281)     birth rejection leaves the node's cached sumWeightedRespon MISSED -
M4   [src/bartcore/moves.hpp:338](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L338)     death rejection never reattaches the orphaned children     caught leafOf matches the derived map after every sweep
M5   [src/bartcore/moves.hpp:569](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L569)     change move's proposal correction sign flipped (both-ordin caught saved dump of forest 0 is mu's, reading the full store
M6   [src/bartcore/moves.hpp:379](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L379)     good-ordinal-rule upper bound off by one (admits a descend caught a wholly stranded tree reports no NaN acceptance probabili
M7   [src/bartcore/moves.hpp:776](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L776)     swap move applies the parent rule to the left child whiche caught crash SIGBUS, no check fired
M8   [src/bartcore/moves.hpp:659](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L659)     categorical gauge check drops the directions == reachable  MISSED -
M9   [src/bartcore/moves.hpp:635](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L635)     ordinal descendant interval off by one on both sides       MISSED -
M10  [src/bartcore/moves.hpp:80](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L80)      branch veto rank never raised above 0                      caught multi-chain rollback leaves valid partitions
C1   [src/bartcore/chain.hpp:1518](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1518)    latent refresh no longer invalidates the vector leaf's U'W MISSED -
C2   [src/bartcore/chain.hpp:1497](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1497)    saved-tree slot ignores savedSlotBase (ring cursor lost)   caught the second run's draws are the three most recent, in order
C3   [src/bartcore/chain.hpp:4207](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4207)    variance forest binds y + meanFits instead of the residual MISSED -
C4   [src/bartcore/chain.hpp:2393](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L2393)    combined variance recompute reads factorByTree with the wr caught variance forest recovers higher s^2 where the truth is noi
C5   [src/bartcore/chain.hpp:651](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L651)     single-forest sigma ladder drops the aft arm               caught aft reduction: uncensored fit is bitwise gaussian
C6   [src/bartcore/chain.hpp:5028](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L5028)    latentScaleAnchor's logistic arm returns probit's constant caught calibration map: the reported scale is the map's product (
C7   [src/bartcore/chain.hpp:4711](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4711)    finalizeTotalFits retires tree 0's fits, not the last tree caught totalFits sums the per-tree fits after a reset sweep
C8   [src/bartcore/chain.hpp:4758](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4758)    subtractTreeFitsFromTotal adds instead of subtracting      caught fuzz [gaussian] seed 1 op 3: totalFits != tree-order sum
C9   [src/bartcore/chain.hpp:3422](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L3422)    column-mask containment gate never refuses                 caught column mask: setState names that refusal the column-mask o
A1   [src/bartcore/sampler.hpp:435](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L435)   recordedDraws forgets the previous run's draws             caught a full store stays full across a partial run
A2   [src/bartcore/sampler.hpp:499](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L499)   savedSlotForDraw drops the oldest-first rotation           caught a partly filled store reads from slot 0
A3   [src/bartcore/sampler.hpp:1733](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L1733)  per-observation session's last-occupant check never fires  caught aggressive per-observation update finalizes
A4   [src/bartcore/sampler.hpp:1749](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L1749)  commit leaves the vacated leaf's occupancy count stale     caught aggressive per-observation update finalizes
A5   [src/bartcore/sampler.hpp:1656](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L1656)  a rolled-back predictor transaction never restores the cod caught rollback restores codes
A6   [src/bartcore/sampler.hpp:765](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L765)   setState ignores the column-mask containment verdict       MISSED -
A7   [src/bartcore/sampler.hpp:785](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L785)   setState restores recordedDraws as a full store            caught saved trees ride along with the state
A8   [src/bartcore/sampler.hpp:575](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L575)   predict writes at the capacity stride, reads at the draw c MISSED -
A9   [src/bartcore/sampler.hpp:277](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L277)   run() pins every chain's saved slot base at 0              caught the second run's draws are the three most recent, in order
D1   [src/bartcore/model.hpp:182](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L182)     constant-leaf marginal's shrinkage term takes the wrong si caught every cell of the exact law is measurable at this draw cou
D2   [src/bartcore/model.hpp:192](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L192)     leaf posterior mean drops the prior precision from its den caught posterior draw mean (actual 0.120015013713005, expected 0.
D3   [src/bartcore/model.hpp:2131](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L2131)    growth probability loses the availability veto             caught single-category node cannot grow
D4   [src/bartcore/model.hpp:2185](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L2185)    ordinal rule prior drops the missing-direction widening    MISSED -
D5   [src/bartcore/model.hpp:2227](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L2227)    categorical group mass is not spread over the emitted cand caught grown categorical root rules follow the exact CGM law at P
D6   [src/bartcore/model.hpp:3560](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L3560)    logistic working response drops the count weight           caught weighted logistic calibrates to the low-quartile weighted 
D7   [src/bartcore/model.hpp:3603](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L3603)    a logistic count swap leaves an inactive row's stale laten MISSED -
D8   [src/bartcore/model.hpp:3314](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L3314)    ordinal category likelihood drops its lower tail           MISSED -
T1   [src/bartcore/tree.hpp:431](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L431)      ancestor cut interval off by one on the left bound         caught crash SIGBUS, no check fired
T2   [src/bartcore/tree.hpp:429](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L429)      splitInterval reads the child direction backwards          caught crash SIGBUS, no check fired
T3   [src/bartcore/tree.hpp:363](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L363)      member-empty leaf demoted to the weight veto rank          caught a member-empty leaf still vetoes with no weights
T4   [src/bartcore/tree.hpp:981](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L981)      death merge drops the right child's weighted response      caught orphan merge sum wz (actual 50, expected 200)
T5   [src/bartcore/tree.hpp:905](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L905)      children's index segments overlap after a partition        caught crash SIGSEGV, no check fired
T6   [src/bartcore/tree.hpp:964](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L964)      birth leaves the right child's leaf statistics stale       caught collapsed subtree takes the plain mean of its leaf paramet
T7   [src/bartcore/tree.hpp:375](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L375)      node depth counted one too deep (prior growth veto shifts) caught root growth probability (actual 0.2375, expected 0.95)
F1   [src/bartcore/facade.hpp:429](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L429)    shape reports the store capacity as the retained draw coun MISSED -
F2   [src/bartcore/facade.hpp:439](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L439)    shape's active-row capability answers the counts probe     caught shape: a gaussian sampler accepts an active-row mask
F3   [src/bartcore/facade.hpp:454](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L454)    setResponse forwards the inverted updateScale flag         MISSED -
F4   [src/bartcore/facade.hpp:514](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L514)    savedSlotForDraw forwarded as the identity                 MISSED -
B1   [src/bartcore/combiner.hpp:894](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L894)  per-forest residual subtracts forest f's own contribution  caught BCF grow-from-root combined fit characteristic value (actu
B2   [src/bartcore/combiner.hpp:896](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L896)  per-forest response is not divided by the amplitude        caught BCF grow-from-root combined fit characteristic value (actu
B3   [src/bartcore/combiner.hpp:897](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L897)  per-forest precision carries m rather than m^2             caught BCF response swap leaves sigma continuous
B4   [src/bartcore/combiner.hpp:944](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L944)  combinedFits accumulates forward instead of last-forest-do caught the combination never absorbs the last forest's product
B5   [src/bartcore/combiner.hpp:911](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L911)  veto precisions lose the near-zero multiplier snap         MISSED -
B6   [src/bartcore/combiner.hpp:1791](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L1791) multinomial forest response subtracts the margin instead o caught K = 2 multinomial reproduces the logistic working-response
B7   [src/bartcore/combiner.hpp:1809](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L1809) multinomial forest precisions ignore the active-row mask   caught an inactive row's composed precision is exactly zero
B8   [src/bartcore/combiner.hpp:1791](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L1791) multinomial PG kappa uses the whole trial count, not half  caught K = 2 multinomial reproduces the logistic working-response
N1   [src/bartcore/chain.hpp:1475](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1475)    accepted move never patches leafOf (map fresh but wrong)   caught crash SIGABRT, no check fired
M3b  [src/bartcore/moves.hpp:280](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L280)     birth rejection restores the node's cached sums as garbage caught dart concentrates splits on signal variables
S2b  [src/bartcore/scan.hpp:159](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L159)      RERUN drop the paired occupancy sentinel on entry[1]       caught crash SIGSEGV, no check fired
S3b  [src/bartcore/scan.hpp:146](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/scan.hpp#L146)      RERUN off-by-one in the prefix scan                        caught crash SIGSEGV, no check fired
G7b  [src/bartcore/grow.hpp:284](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L284)      RERUN cut buffer sized numCuts, not 2*numCuts              caught crash SIGTRAP, no check fired
M7b  [src/bartcore/moves.hpp:776](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L776)     RERUN swap applies the parent rule to the left child regar caught joint per-observation update finalizes
F5   [src/bartcore/facade.hpp:418](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L418)    shape reports the chain count as the forest count          caught shape: constant gaussian is single-forest
F6   [src/bartcore/facade.hpp:619](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L619)    per-forest weights installed on forest 0 whatever was aske MISSED -
F7   [src/bartcore/facade.hpp:519](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L519)    saved-tree reads take the reported forest's slot regardles MISSED -

## Ladder audit
Probe: a 7th ResponseFamily enumerator in a staged copy; tests/cpp rebuilt and R_interface_bartcore.cpp
compiled -fsyntax-only, -Wall -Wextra.  Exactly THREE sites warned (-Wswitch): [src/bartcore/chain.hpp:582](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L582),
R_interface:2268, [src/bartcore/chain.hpp:5014](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L5014).  Four `default:` ladders stayed silent, plus one open-coded chain no -Wswitch reaches
(only two of the four are in src/bartcore; the memo's list spans both files):
  [src/bartcore/chain.hpp:756](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L756)  AmplitudeSpec ctor        aft/ord/nbinom/new  -> GaussianResponse
  [src/bartcore/chain.hpp:5026](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L5026) latentScaleAnchor         gauss/aft/ord/nb/new -> scaledResponseSd()
  R_interface:2291 defaultNodeScale        gaussian/aft/new    -> 0.5 (gaussian's)
  R_interface:2812 validateResponseSupport gaussian/aft/new    -> NO validation
  R_interface:6254 computeWorkingResponse  if/else, no switch  -> latent[i]
[src/bartcore/chain.hpp:756](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L756) and [src/bartcore/chain.hpp:5026](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L5026) take spec.family from the K-forest route, which [src/R_interface_bartcore.cpp:2268](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L2268) refuses for aft/ordinal/nbinom first, so
those arms are gaussian-only today.  [src/R_interface_bartcore.cpp:2812](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L2812) is the worst door: its own comment enumerates the harm it prevents
(a negative nbinom count underflowing into a ~1.8e19 allocation, "an uncatchable crash"). The de-facto 7th
family already exists - [src/bartcore/chain.hpp:854](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L854) sets `family_ = logistic` for multinomial commenting "family() is not
read on this path", yet SamplerShape::family ([src/bartcore/facade.hpp:437](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L437)) reports it and the bridge branches on it at
six sites (unresolved: [src/R_interface_bartcore.cpp:2729](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L2729), unresolved: [src/R_interface_bartcore.cpp:4611](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L4611), unresolved: [src/R_interface_bartcore.cpp:4662](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L4662), unresolved: [src/R_interface_bartcore.cpp:4702](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L4702), unresolved: [src/R_interface_bartcore.cpp:4877](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L4877), unresolved: [src/R_interface_bartcore.cpp:4949](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/R_interface_bartcore.cpp#L4949) - the exact line within the bridge that each reference names could not be re-derived), each unreachable only because a SEPARATE multi-forest refusal stands in
front.  VERDICT the doors are shut by refusals in another file, not by the type system, and one family already
travels under another's name.

## Assertion audit
R compiles packages with -DNDEBUG (Makeconf:174), so every runtime assert() in src/bartcore is DEAD in the
shipped library; tests/cpp never defines it, so the suite runs them LIVE.  22 static_asserts are compile-time
- no gap.  The six runtime asserts, with their release-build counterpart:
  [src/bartcore/chain.hpp:4425](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4425) recoveredFactors || bottomNodesAreOccupied - YES: recoverVarianceLeafValuesBelow's own
    numObservations() > 0 test returns the 1.0 identity instead of dereferencing.
  [src/bartcore/chain.hpp:4494](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4494)/[src/bartcore/chain.hpp:4522](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4522) leaf[j] < muByTree[t].size() - NO; its own comment concedes a stale index reads inside
    mu's CAPACITY, returning a stale value rather than faulting, invisible to ASAN too.
  [src/bartcore/chain.hpp:4628](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4628) leaf[j] < numNodes && isBottom - PARTIAL: the fused pass declines a map FLAGGED stale
    ([src/bartcore/chain.hpp:4617](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4617)), not one that is fresh and WRONG.
  [src/bartcore/grow.hpp:249](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L249), [src/bartcore/model.hpp:2224](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L2224) numPresent in [2, numReachable] - NO, but the failure is a finite wrong prior
    mass, not a fault.
EMPIRICAL: N1 makes leafOf fresh-but-wrong.  Plain
build SIGABRTs in the GROW suite, before any leafOf test runs; rebuilt -DNDEBUG (R's own flag) it is rc 1 with
53 failures, first "leafOf matches the derived map after every sweep" ([tests/cpp/test_moves.cpp:1456](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/tests/cpp/test_moves.cpp#L1456)); M10, the only
other SIGABRT, also reddens (70). VERDICT no catch here depends on an assertion the release build lacks.  The
residual exposure is that the SHIPPED library has no runtime guard on the (mu, leafOf) pairing at all - only a
test does.  The sentinel shape, milder.

## Findings
1 BLOCKER [src/bartcore/chain.hpp:1518](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1518) (also [src/bartcore/chain.hpp:1578](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1578),[src/bartcore/chain.hpp:1616](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1616),[src/bartcore/chain.hpp:1665](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1665),[src/bartcore/chain.hpp:1678](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L1678)) [C1].  Deleting the invalidateStatistics the latent
  refresh owes the vector leaf's U'WU cache passes.  Every LinearGaussianLeaf/GPGaussianLeaf fixture in
  tests/cpp is gaussian, whose workingWeightsVaryPerSweep() is false, so the branch all five call sites live
  in is NEVER ENTERED; the primitive is unit-tested (test_model:1393), the wiring is not.  Reachable and
  unguarded elsewhere: [src/bartcore/facade.hpp:816](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L816) has no family gate, R/spec.R refuses linear/GP leaves only on the
  multinomial and K-forest routes, and no tinytest pairs a linear/gp node prior with a non-gaussian family.
  ASSERTION one end-to-end fixture per vector leaf model under a workingWeightsVaryPerSweep family, its draws
  compared with a cache-disabled run.  FIX agent-fix (test); refusing the cross at the bridge instead is
  VD-judgement.
2 BLOCKER facade.hpp forwarding [F1,F3,F4,F6,F7].  5 of 7 facade mutations survive: numSavedDraws reporting
  capacity ([src/bartcore/facade.hpp:429](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L429)), setResponse forwarding !updateScale ([src/bartcore/facade.hpp:454](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L454)), savedSlotForDraw as the identity ([src/bartcore/facade.hpp:514](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L514)),
  setForestWeights on forest 0 regardless ([src/bartcore/facade.hpp:619](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L619)), savedTree on forest 0 regardless ([src/bartcore/facade.hpp:519](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/facade.hpp#L519)); only the two
  shape() fields were caught.  Every other test drives Sampler<L> DIRECTLY, so the 36-virtual facade - the one
  dispatch layer between the shipped flat C API and the engine - is exercised only through shape(), while the
  ENGINE twins of three of these WERE caught (A2, A7). ASSERTION a facade conformance test: for each virtual
  with a selecting argument (forestIndex, slot, chainNum, updateScale), one call through SamplerBase whose
  answer differs per argument value.  FIX agent-fix (test).
3 MAJOR [src/bartcore/chain.hpp:4207](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L4207) [C3].  Binding y + meanFits for y - meanFits - the variance forest fits the wrong
  quantity - passes: the only VALUE-level variance assertion (test_model:689-725) uses y = s(x)*N(0,1) with NO
  MEAN FUNCTION, so y-f and y+f differ by 2f ~ 0, and its bound is loose. MAJOR not BLOCKER - tinytest
  backstops it ([inst/tinytest/test-heteroscedastic.R:11-42](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/inst/tinytest/test-heteroscedastic.R#L11-L42)) - but [docs/design/heteroscedastic.md:465](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/docs/design/heteroscedastic.md#L465) states gate (d) as "recovers f(x)
  AND s(x)" and the C++ fixture dropped f(x).  ASSERTION the same test on a strong non-constant mean, s^2(x)
  tracking the noise and NOT |f(x)|.  FIX agent-fix (test).
4 MAJOR [src/bartcore/sampler.hpp:575](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L575) (identically [src/bartcore/sampler.hpp:605](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L605), [src/bartcore/sampler.hpp:679](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L679)) [A8].  Predict's destination taking the CAPACITY stride while
  the loop runs over filledSavedDraws() passes; regressed that is an out-of-bounds heap WRITE for a partially
  filled store and any chain c >= 1.  Every saved-tree test is single-chain (test_state:255 leaves numChains
  at 1, so c*numDraws == c*capacity == 0), so the sanitizer cannot help either.  ASSERTION repeat the
  partial-fill section with numChains >= 2, checking chain 1's slab starts at numDraws*slab, for all three
  readers.  FIX agent-fix (test).
5 MAJOR [src/bartcore/moves.hpp:635](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L635) and [src/bartcore/moves.hpp:659](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L659) [M9, M8].  ordinalRuleIsValid's descendant interval off by one on both sides,
  and categoricalSubtreeIsValid's gauge dropping `directions == reachable`, both survive.  ruleIsValid,
  ordinalRuleIsValid, categoricalSubtreeIsValid(Wide) and findGoodOrdinalRules have ZERO direct callers in
  tests/cpp; reached only as a side condition of swap/change, a wrong verdict moves a proposal between "no-op"
  and "scored", which no structural check sees.  ASSERTION direct unit tests over hand-built trees - at the
  ancestor bound accepted and past it refused; a mask equal to the reachable set refused, a strict nonempty
  subset accepted.  FIX agent-fix.
6 MAJOR [src/bartcore/grow.hpp:222](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/grow.hpp#L222) [G5].  -log(numPredictors) for -log(numAvailable) - identical at a root, wrong at every
  node below - survives: all three grow-from-root LAW tests measure ROOT rules only, and below the root the
  suite asserts legality and gauge, never a probability.  ASSERTION one conditional law test at depth >= 1
  (fix the root rule, chi-square a child's realized rules against the exact law).  FIX agent-fix (test).
7 MAJOR [src/bartcore/model.hpp:3314](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L3314) [D8].  Dropping the lower tail from OrdinalResponse::computeLogLikelihood survives -
  the EXPORTED per-observation log-likelihood channel of every ordinal fit, the shape the 08-24 value scan
  already found four times.  ASSERTION pin it against an independently coded Phi difference on a small ordinal
  fixture, as the file already does for the cutpoint acceptance. FIX agent-fix (test).
8 MAJOR [src/bartcore/model.hpp:2185](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L2185) [D4].  Dropping the `+ log 2` that ruleForVariableLogProbability owes a missing-bearing
  ordinal column survives. It enters treeLogProbability and both MH ratios but cancels on a same-variable
  redraw - which is why the grow-side twin (G1) is caught and this is not: nothing measures a ratio moving
  BETWEEN a missing-bearing and a plain column.  ASSERTION pin ruleForVariableLogProbability against a
  hand-computed value on both column kinds.  FIX agent-fix (test).
9 MAJOR [src/bartcore/model.hpp:3603](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/model.hpp#L3603) [D7].  Deleting LogisticResponse::setWeights' cold start of INACTIVE rows survives -
  the line that IS the logistic weight channel's landing claim; no landed test constructs a count swap while a
  mask is installed.  ASSERTION install a mask, swap the counts, clear the mask, check the reactivated rows'
  latents equal a cold start against the NEW counts. FIX agent-fix (test).
10 MINOR [src/bartcore/combiner.hpp:911](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L911) [B5].  formForestVetoWeights losing the near-zero multiplier snap survives;
  testBCFZeroMultiplierSnap covers the response half ([src/bartcore/combiner.hpp:894-897](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/combiner.hpp#L894-L897), all three caught), not the veto half the
  empty-leaf veto reads. ASSERTION extend it: a snapped row's veto precision is exactly 0.0.  FIX agent-fix
  (test).
11 MINOR diagnosis, not detection.  16 of 63 catches were crashes naming nothing (stdout lost to buffering),
  so CI reddens with no FAIL line. ASSERTION none; `setvbuf(stdout, NULL, _IOLBF, 0)` in tests/cpp/main.cpp
  makes a crashing run name the check it died in.  FIX agent-fix (one line).

## Real defects found incidentally
None user-visible in pristine code; two code-health items, with reproductions.
R1 [src/bartcore/moves.hpp:280-281](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/moves.hpp#L280-L281) - the birth-rejection restore of the parent's sumWeights/sumWeightedResponse is DEAD:
  Tree::birth ([src/bartcore/tree.hpp:953](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/tree.hpp#L953)) never writes the parent's sums.  Reproduction: M3 deletes one line and the suite
  passes (equivalent mutant); M3b replaces both with 1.0e6 and the suite reddens, so the fields are read later
  but the restore changes nothing, while the comment above reads as load-bearing.  FIX agent-fix (delete, fix
  comment).
R2 [src/bartcore/sampler.hpp:765](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/sampler.hpp#L765) - `bool allValid = columnMaskOk;` is redundant with Chain::stateIsValid's own per-forest
  ([src/bartcore/chain.hpp:3229](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L3229)) and per-variance-tree ([src/bartcore/chain.hpp:3318](https://github.com/vdorie/dbarts/blob/9cebb35221ff0d932f126c1a8f710eb464fbc608/src/bartcore/chain.hpp#L3318)) columnMaskSubtreeIsValid, and *columnMaskRefused reads
  columnMaskOk directly, so `true` is behaviourally equivalent (A6, missed).  FIX defer.

## Sanitizer leg
Pristine staged tree, `make clean && make OPT="-O2 -g -fsanitize=address,undefined"`, run once under
ASAN_OPTIONS=detect_container_overflow=0: build rc 0, run rc 0, "all tests passed", ZERO "runtime error" or
"AddressSanitizer" lines.  Clean.
