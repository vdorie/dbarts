// A C++ consumer of the flat C API whose sweep hook and draw callback throw,
// against the header's rule, compiled by test-capi.R with R CMD SHLIB against
// the installed headers: dbarts' safety net catches an exception thrown in this
// library and hands it back to R as an error carrying its what(), from the
// calling thread and from a worker. Both libraries share one C++ runtime here,
// which the net needs and the header does not promise.

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#define DBARTS_USE_STUBS
#include <dbarts/dbarts.h>

#include <stdexcept>

namespace {

// the last chain throws at its third sweep, or at its second draw
int throwingSweep(void* context, dbarts_chain_view* view, size_t sweepIndex,
                  int) {
  size_t last = *static_cast<size_t*>(context) - 1;
  if (view->chainIndex == last && sweepIndex == 2)
    throw std::runtime_error("the sweep hook threw");
  return 0;
}

int throwingDraw(void* context, const dbarts_draw* draw) {
  size_t last = *static_cast<size_t*>(context) - 1;
  if (draw->chainIndex == last && draw->drawIndex == 1)
    throw std::runtime_error("the draw callback threw");
  return 0;
}

size_t numChains = 0;

}  // namespace

// which: 1 registers the throwing sweep hook, 2 the throwing draw callback,
// 0 clears both
extern "C" SEXP throwing_register(SEXP ptrExpr, SEXP whichExpr) {
  dbarts_sampler* sampler =
    static_cast<dbarts_sampler*>(R_ExternalPtrAddr(ptrExpr));
  int which = Rf_asInteger(whichExpr);
  numChains = dbarts_sampler_numChains(sampler);
  dbarts_sampler_setSweepHook(sampler, which == 1 ? &throwingSweep : nullptr,
                              &numChains);
  dbarts_sampler_setDrawCallback(
    sampler, which == 2 ? &throwingDraw : nullptr, &numChains, DBARTS_DRAW_ALL);
  return R_NilValue;
}

extern "C" SEXP throwing_run(SEXP ptrExpr, SEXP numSamplesExpr) {
  dbarts_sampler_run(static_cast<dbarts_sampler*>(R_ExternalPtrAddr(ptrExpr)),
                     0, static_cast<size_t>(Rf_asInteger(numSamplesExpr)),
                     nullptr);
  return R_NilValue;
}
