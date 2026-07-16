# capi-dispatch-table

status: investigation in flight (queued 2026-07-15 by VD; the design memo
lands in this file)

VD's ask: do a design pass on the C ABI; see if it is worth breaking with
R's mechanism for passing C function pointers by symbol
(R_RegisterCCallable/R_GetCCallable per function) and instead doing
something like a single function that returns a function table, indexed
by a version string.

Motivating failure mode: the getTrees ABI break (2e2b1c9, restored
a73ca50). Consumers declare the signatures themselves, so a signature
change under an unchanged symbol name is a silent stack-garbage break; a
versioned table can make the mismatch detectable at load time.

Priors: docs/plans/c-api-growth.md (the structSize results-struct growth
mechanism and by-name state blocks), docs/design/public-surface.md
section 6 (the stan4bart usage audit). window: pre-release - the
dispatch mechanism is only cheap to change before 1.0-0 freezes the
dbarts.h compatibility contract.
