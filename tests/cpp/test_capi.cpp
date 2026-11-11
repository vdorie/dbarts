// The flat C API's per-draw and per-sweep surface: the shipped dbarts_draw
// against the engine's DrawInfo, the registration semantics
// dbarts_sampler_setDrawCallback and dbarts_sampler_setSweepHook have, and the
// chain view's operations. The entry point itself speaks SEXP and is linked
// only by the package build, so what runs here is the state and the adapter it
// delegates to - one line apiece there, and all of the behavior.

#include "common.hpp"

#include "R_interface_bartcore_common.hpp"

using bartcore_bridge::fillShippedDraw;
using bartcore_bridge::ShippedDrawHook;
using bartcore_bridge::ShippedSweepHook;

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

// Every channel the engine fills reaches the callback, at the library's own
// structSize: a field the adapter forgets keeps the value the stack copy found.
void testAdapterCopiesEveryField() {
  double values[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
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
  info.latents = values + 8;
  info.sigma = 0.25;
  info.k = 2.5;
  info.shape = std::numeric_limits<double>::quiet_NaN();
  info.residualDf = 4.5;

  dbarts_draw draw;
  std::memset(&draw, 0xa5, sizeof draw);
  fillShippedDraw(draw, info);

  bool copied = draw.structSize == sizeof(dbarts_draw) &&
    draw.chainIndex == 3 && draw.drawIndex == 5 && draw.numObservations == 7 &&
    draw.numTestObservations == 9 && draw.numPredictors == 11 &&
    draw.numFittedValuesPerObservation == 2 &&
    draw.numVariableCountForests == 4 &&
    draw.numForests == 6 && draw.numAmplitudes == 8 &&
    draw.numOrdinalThresholds == 10 &&
    draw.train == values + 0 && draw.test == values + 1 &&
    draw.varianceFits == values + 2 && draw.varianceTestFits == values + 3 &&
    draw.forestFits == values + 4 && draw.glue == values + 5 &&
    draw.splitProbabilities == values + 6 && draw.logLikelihood == values + 7 &&
    draw.ordinalThresholds == nullptr && draw.varcount == counts &&
    draw.latents == values + 8 &&
    draw.sigma == 0.25 && draw.k == 2.5 && std::isnan(draw.shape) &&
    draw.residualDf == 4.5;
  check(copied, "capi draw: the adapter copies every channel and scalar");
}

// register, clear, re-register: the whole of the setter.
void testRegistration() {
  ShippedDrawHook hook;
  check(hook.engineHook().fn == nullptr,
        "capi draw: nothing is registered by default");

  CallLog first;
  hook.set(&logCall, &first, DBARTS_DRAW_TRAIN | DBARTS_DRAW_LATENTS);
  DrawHook engine = hook.engineHook();
  check(engine.fn != nullptr &&
          engine.channels == (bartcore::drawTrain | bartcore::drawLatents),
        "capi draw: a registration installs a hook for the channels named");
  DrawInfo info;
  info.chainIndex = 2;
  check(engine.fn(engine.context, &info) == 0 && first.calls == 1 &&
          first.context == &first && first.chainIndex == 2 &&
          first.structSize == sizeof(dbarts_draw),
        "capi draw: the hook calls the registered function with its context");

  // a second registration REPLACES both halves
  CallLog second;
  second.returnValue = 1;
  hook.set(&logCall, &second, DBARTS_DRAW_ALL);
  engine = hook.engineHook();
  check(engine.fn(engine.context, &info) == 1 && second.calls == 1 &&
          first.calls == 1 && engine.channels == bartcore::drawAllChannels,
        "capi draw: a second registration replaces the function, context and "
        "channels");

  // a null function clears, and drops the context with it
  hook.set(nullptr, &second, DBARTS_DRAW_ALL);
  check(hook.fn == nullptr && hook.context == nullptr && hook.channels == 0 &&
          hook.engineHook().fn == nullptr,
        "capi draw: a null function clears the hook, its context and channels");
  // and a cleared hook installs nothing, so refuseCall is unreachable
  hook.set(&refuseCall, nullptr, DBARTS_DRAW_ALL);
  hook.set(nullptr, nullptr, 0);
  DrawHook cleared = hook.engineHook();
  check(cleared.fn == nullptr && cleared.context == nullptr,
        "capi draw: a cleared hook is empty on both halves");
}

/// The engine's side of a chain, recording what a view's operations hand it.
struct ChainSpy final : ChainView {
  const double* offset = nullptr;
  int offsetWrites = 0, sigmaWrites = 0;
  double sigma = 0.0;
  void setOffsetHoldingScale(const double* values) override {
    offset = values;
    ++offsetWrites;
  }
  void setSigma(double value) override {
    sigma = value;
    ++sigmaWrites;
  }
};

/// What a probing hook saw, and each operation's status in call order.
struct ViewProbe {
  dbarts_chain_view seen;
  std::size_t sweepIndex = 0;
  int isBurnIn = -1;
  std::vector<int> statuses;
};

int probeSweep(void* context, dbarts_chain_view* view, std::size_t sweepIndex,
               int isBurnIn) {
  ViewProbe& probe = *static_cast<ViewProbe*>(context);
  probe.seen = *view;
  probe.sweepIndex = sweepIndex;
  probe.isBurnIn = isBurnIn;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const double bad[4] = { 1.0, nan, 3.0, 4.0 };
  const double badLast[4] = { 1.0, 2.0, 3.0, nan };
  const double good[4] = { 1.0, 2.0, 3.0, 4.0 };
  probe.statuses = { view->setOffset(view, bad), view->setOffset(view, badLast),
                     view->setOffset(view, good), view->setSigma(view, 0.0),
                     view->setSigma(view, inf), view->setSigma(view, 2.0) };
  return 7;
}

// The view a registered hook is handed: its chain, its own buffer, and each
// operation's status, a refused value reaching nothing.
void testSweepHookView() {
  ShippedSweepHook hook;
  check(hook.prepareRun(3, 4, false).fn == nullptr,
        "capi sweep: nothing is registered by default");

  ViewProbe probe;
  hook.set(&probeSweep, &probe, 3, 4);
  ChainSpy chain;
  SweepHook engine = hook.prepareRun(3, 4, false);
  int result = engine.fn(engine.context, chain, 2, 7, true);
  const double* own = hook.chainOffsets[2].data();
  check(result == 7 && probe.seen.structSize == sizeof(dbarts_chain_view) &&
          probe.seen.chainIndex == 2 && probe.seen.numObservations == 4 &&
          probe.sweepIndex == 7 && probe.isBurnIn == 1,
        "capi sweep: the hook sees its chain, sweep and burn-in flag, and its "
        "return goes back to the engine");
  check(probe.statuses == std::vector<int>{ DBARTS_INVALID_ARGUMENT,
                                            DBARTS_INVALID_ARGUMENT,
                                            DBARTS_OK,
                                            DBARTS_INVALID_ARGUMENT,
                                            DBARTS_INVALID_ARGUMENT,
                                            DBARTS_OK } &&
          chain.offsetWrites == 1 && chain.offset == own && own[1] == 2.0 &&
          own[3] == 4.0 && chain.sigmaWrites == 1 && chain.sigma == 2.0,
        "capi sweep: a non-finite offset, in any row, or a sigma not finite "
        "and positive is refused untouched; a write lands in the chain's own "
        "buffer");

  ChainSpy pinned;
  engine = hook.prepareRun(3, 4, true);
  engine.fn(engine.context, pinned, 0, 0, false);
  check(probe.statuses.back() == DBARTS_UNSUPPORTED &&
          pinned.sigmaWrites == 0 && probe.isBurnIn == 0 &&
          pinned.offset == hook.chainOffsets[0].data(),
        "capi sweep: a pinned sigma answers unsupported and writes nothing");

  hook.set(nullptr, &probe, 3, 4);
  check(hook.fn == nullptr && hook.context == nullptr &&
          hook.chainOffsets.empty() &&
          hook.prepareRun(3, 4, false).fn == nullptr,
        "capi sweep: a null function clears the hook, its context and the "
        "buffers");
}

}  // namespace

void runCapiTests() {
  testAdapterCopiesEveryField();
  testRegistration();
  testSweepHookView();
  printf("ok: flat C API per-draw struct, adapter and registration\n");
  printf("ok: flat C API sweep hook registration and chain view\n");
}
