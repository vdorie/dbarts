// The flat C API's per-draw surface: the shipped dbarts_draw against the
// engine's DrawInfo, and the registration semantics dbarts_sampler_setDrawCallback
// has. The entry point itself speaks SEXP and is linked only by the package
// build, so what runs here is the state and the adapter it delegates to - one
// line apiece there, and all of the behaviour.

#include "common.hpp"

#include "R_interface_bartcore_common.hpp"

using bartcore_bridge::fillShippedDraw;
using bartcore_bridge::ShippedDrawHook;

namespace {

struct CallLog {
  int calls = 0;
  void* context = nullptr;
  std::size_t structSize = 0;
  std::size_t chainIndex = 0;
  int returnValue = 0;
};

int logCall(void* context, const dbarts_draw* draw) {
  CallLog& log = *static_cast<CallLog*>(context);
  ++log.calls;
  log.context = context;
  log.structSize = draw->structSize;
  log.chainIndex = draw->chainIndex;
  return log.returnValue;
}

int refuseCall(void*, const dbarts_draw*) {
  check(false, "capi draw: a cleared hook must not be called");
  return 0;
}

// The shipped struct mirrors the engine's field order. Both are standard
// layout, so the two offset runs are directly comparable: a field moved on one
// side alone shows up as an inversion here, which no compile-time assert on
// either struct by itself would see.
void testFieldOrderMirrorsEngine() {
  struct Pair { std::size_t shipped, engine; };
  const Pair pairs[] = {
#define PAIR(field) { offsetof(dbarts_draw, field), offsetof(DrawInfo, field) }
    PAIR(chainIndex), PAIR(drawIndex), PAIR(numObservations),
    PAIR(numTestObservations), PAIR(numPredictors), PAIR(numReportedLocations),
    PAIR(numVariableCountForests), PAIR(numForests), PAIR(numAmplitudes),
    PAIR(numOrdinalThresholds), PAIR(train), PAIR(test), PAIR(varianceFits),
    PAIR(varianceTestFits), PAIR(forestFits), PAIR(glue),
    PAIR(splitProbabilities), PAIR(logLikelihood), PAIR(ordinalThresholds),
    PAIR(varcount), PAIR(sigma), PAIR(k), PAIR(dispersion), PAIR(residualDf)
#undef PAIR
  };
  bool ordered = true;
  for (std::size_t i = 1; i < std::size(pairs); ++i)
    ordered &= pairs[i].shipped > pairs[i - 1].shipped &&
               pairs[i].engine > pairs[i - 1].engine;
  check(ordered, "capi draw: shipped field order mirrors the engine struct's");
  // structSize leads, and it is the only field the engine struct does not have
  check(offsetof(dbarts_draw, structSize) == 0 &&
          pairs[0].shipped == sizeof(std::size_t),
        "capi draw: structSize leads the shipped struct alone");
}

// Every channel the engine fills reaches the callback, at the library's own
// structSize: a field the adapter forgets keeps the value the stack copy found.
void testAdapterCopiesEveryField() {
  double values[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  std::uint32_t counts[2] = { 11, 12 };
  DrawInfo info;
  info.chainIndex = 3;
  info.drawIndex = 5;
  info.numObservations = 7;
  info.numTestObservations = 9;
  info.numPredictors = 11;
  info.numReportedLocations = 2;
  info.numVariableCountForests = 4;
  info.numForests = 6;
  info.numAmplitudes = 8;
  info.numOrdinalThresholds = 10;
  info.train = values + 0;
  info.test = values + 1;
  info.varianceFits = values + 2;
  info.varianceTestFits = values + 3;
  info.forestFits = values + 4;
  info.glue = values + 5;
  info.splitProbabilities = values + 6;
  info.logLikelihood = values + 7;
  info.ordinalThresholds = nullptr;  // an absent channel stays null
  info.varcount = counts;
  info.sigma = 0.25;
  info.k = 2.5;
  info.dispersion = std::numeric_limits<double>::quiet_NaN();
  info.residualDf = 4.5;

  dbarts_draw draw;
  std::memset(&draw, 0xa5, sizeof draw);
  fillShippedDraw(draw, info);

  bool copied = draw.structSize == sizeof(dbarts_draw) &&
    draw.chainIndex == 3 && draw.drawIndex == 5 && draw.numObservations == 7 &&
    draw.numTestObservations == 9 && draw.numPredictors == 11 &&
    draw.numReportedLocations == 2 && draw.numVariableCountForests == 4 &&
    draw.numForests == 6 && draw.numAmplitudes == 8 &&
    draw.numOrdinalThresholds == 10 &&
    draw.train == values + 0 && draw.test == values + 1 &&
    draw.varianceFits == values + 2 && draw.varianceTestFits == values + 3 &&
    draw.forestFits == values + 4 && draw.glue == values + 5 &&
    draw.splitProbabilities == values + 6 && draw.logLikelihood == values + 7 &&
    draw.ordinalThresholds == nullptr && draw.varcount == counts &&
    draw.sigma == 0.25 && draw.k == 2.5 && std::isnan(draw.dispersion) &&
    draw.residualDf == 4.5;
  check(copied, "capi draw: the adapter copies every channel and scalar");
}

// register, clear, re-register: the whole of the setter.
void testRegistration() {
  ShippedDrawHook hook;
  check(hook.engineHook().fn == nullptr,
        "capi draw: nothing is registered by default");

  CallLog first;
  hook.set(&logCall, &first);
  DrawHook engine = hook.engineHook();
  check(engine.fn != nullptr, "capi draw: a registration installs a hook");
  DrawInfo info;
  info.chainIndex = 2;
  check(engine.fn(engine.context, &info) == 0 && first.calls == 1 &&
          first.context == &first && first.chainIndex == 2 &&
          first.structSize == sizeof(dbarts_draw),
        "capi draw: the hook calls the registered function with its context");

  // a second registration REPLACES both halves
  CallLog second;
  second.returnValue = 1;
  hook.set(&logCall, &second);
  engine = hook.engineHook();
  check(engine.fn(engine.context, &info) == 1 && second.calls == 1 &&
          first.calls == 1,
        "capi draw: a second registration replaces the function and context");

  // a null function clears, and drops the context with it
  hook.set(nullptr, &second);
  check(hook.fn == nullptr && hook.context == nullptr &&
          hook.engineHook().fn == nullptr,
        "capi draw: a null function clears the hook and its context");
  // and a cleared hook installs nothing, so refuseCall is unreachable
  hook.set(&refuseCall, nullptr);
  hook.set(nullptr, nullptr);
  DrawHook cleared = hook.engineHook();
  check(cleared.fn == nullptr && cleared.context == nullptr,
        "capi draw: a cleared hook is empty on both halves");
}

}  // namespace

void runCapiTests() {
  testFieldOrderMirrorsEngine();
  testAdapterCopiesEveryField();
  testRegistration();
  printf("ok: flat C API per-draw struct, adapter and registration\n");
}
