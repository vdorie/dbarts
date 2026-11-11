#ifndef R_INTERFACE_EXCEPTIONS_HPP
#define R_INTERFACE_EXCEPTIONS_HPP

// How a C++ exception becomes an R error at the package's entrances. Kept
// apart from the bridge internals so the .Call registration table can apply
// it without compiling the engine.

#include <cstddef>     // size_t
#include <cstdio>      // snprintf
#include <exception>   // exception
#include <type_traits> // is_void_v
#include <utility>     // forward

#include <external/Rinternals.h> // SEXP, R_ContinueUnwind, Rf_error

namespace bartcore_bridge {

/// R's continuation for an error longjmp caught at the LEAF that raised it -
/// today the one call out to a consumer's code, a registered draw callback.
/// R_UnwindProtect's cleanup throws this instead of letting R resume the jump,
/// so every C++ frame between that leaf and the entry point unwinds with its
/// destructors, and the entry point resumes the jump with R_ContinueUnwind.
///
/// A longjmp cannot be undone, which is why the protection sits at the leaf
/// rather than at the entry: an entry-level R_UnwindProtect would catch the
/// jump only after the frames beneath it were already abandoned.
struct UnwindJump {
  SEXP continuation;
};

/// The process's one unwind continuation, created and preserved on first use.
/// R fills it on the jump and R_ContinueUnwind consumes it immediately, and
/// only the thread that entered a run ever installs the protection, so one
/// token serves every entrance without per-call allocation or protection-stack
/// bookkeeping.
inline SEXP unwindContinuation() {
  static SEXP token = []() -> SEXP {
    SEXP made = R_MakeUnwindCont();
    R_PreserveObject(made);
    return made;
  }();
  return token;
}

/// A failed call's message, copied out of the exception that carried it. Sized
/// for the engine's own sentences; a longer one truncates rather than
/// allocating, since this path is already failing.
struct CapturedError {
  char message[512];
  bool failed = false;
};

/// Runs \p body, recording a C++ exception escaping it in \p into instead of
/// letting it out. An exception must never cross into R's own C frames - they
/// run no destructors and keep context bookkeeping an unwind past them
/// abandons, and nothing above catches, so the session aborts.
///
/// The message is copied out and the handler LEFT before any caller raises: a
/// longjmp out of a live catch block would strand the exception on the
/// thread's caught-exception stack, so the raise happens where nothing is in
/// flight and every owner the unwind reached is already gone. An UnwindJump is
/// NOT an error and passes straight through: it belongs to the entry point,
/// which resumes the jump it carries.
template <typename Body>
void captureExceptions(CapturedError& into, Body&& body) {
  try {
    body();
  } catch (const UnwindJump&) {
    throw;
  } catch (const std::exception& error) {
    std::snprintf(into.message, sizeof(into.message), "%s", error.what());
    into.failed = true;
  } catch (...) {
    std::snprintf(into.message, sizeof(into.message), "unknown C++ exception");
    into.failed = true;
  }
}

/// captureExceptions, then the raise, labeled with \p caller. For the flat C
/// API, whose errors are read by a compiled consumer that needs to know which
/// entry point refused.
template <typename Body>
void callConvertingExceptions(const char* caller, Body&& body) {
  CapturedError error;
  captureExceptions(error, std::forward<Body>(body));
  if (error.failed) Rf_error("%s: %s", caller, error.message);
}

/// Runs \p body and returns its value, with every C++ outcome turned into the
/// R one: an escaping exception becomes an R error carrying its message alone,
/// and an UnwindJump resumes the jump it carries. Both are raised once the
/// handler is left and every frame beneath has unwound, so the owners in those
/// frames are freed and nothing is in flight at the longjmp.
///
/// This is the one conversion on the R side. Every .Call entry is registered
/// through GuardedEntry, which runs the entry under it, and unwindProtect runs
/// its closure under it, because R_UnwindProtect is itself a C frame. A
/// caller that must act between the failure and the raise (free a buffer by
/// hand, hand back R's stream, release a sampler left part-way, add context to
/// the message) still captures for itself and raises after.
template <typename Body>
SEXP callRaisingExceptions(Body&& body) {
  SEXP result = R_NilValue;
  SEXP continuation = NULL;
  CapturedError error;
  try {
    captureExceptions(error, [&]() { result = body(); });
  } catch (const UnwindJump& jump) {
    continuation = jump.continuation;
  }
  if (continuation != NULL) R_ContinueUnwind(continuation); // does not return
  if (error.failed) Rf_error("%s", error.message);
  return result;
}

/// The registration form of a .Call entry \p Entry: the same signature, the
/// body run under callRaisingExceptions. R calls an entry from its own C
/// frames, so an exception that reached them would abort the session; wrapping
/// at registration rather than in each body is what makes "no exception
/// leaves a .Call entry" hold for every entry, the ones added later included.
/// Costs one handler per call, which takes nothing on the path that does not
/// throw.
template <auto Entry>
struct GuardedEntry;

template <typename... Args, SEXP (*Entry)(Args...)>
struct GuardedEntry<Entry> {
  static SEXP call(Args... args) {
    return callRaisingExceptions([&]() -> SEXP { return Entry(args...); });
  }
};

/// An entry point's name as a template argument, so the flat C API's guard
/// can label its error without the entry passing one.
template <std::size_t N>
struct EntryName {
  char text[N];
  constexpr EntryName(const char (&name)[N]) {
    for (std::size_t i = 0; i < N; ++i) text[i] = name[i];
  }
};

/// The registration form of a flat C API entry \p Entry, named \p Name: the
/// same signature, the body run under callConvertingExceptions, so an
/// exception escaping any entry becomes an R error naming it and carrying its
/// what(), and an UnwindJump resumes the jump it carries. A compiled
/// consumer calls an entry from its own frames, which no exception may cross.
/// A value-returning entry that raises returns nothing, the raise being a
/// jump; this is where an engine reporting errors without R would return a
/// negative status and keep the message instead.
template <auto Entry, EntryName Name, typename Signature = decltype(Entry)>
struct GuardedApiEntry;

template <auto Entry, EntryName Name, typename Result, typename... Args>
struct GuardedApiEntry<Entry, Name, Result (*)(Args...)> {
  static Result call(Args... args) {
    if constexpr (std::is_void_v<Result>) {
      run([&]() { Entry(args...); });
    } else {
      Result result{};
      run([&]() { result = Entry(args...); });
      return result;
    }
  }

private:
  /// callConvertingExceptions, then the jump an UnwindJump carries, resumed
  /// once its handler is left.
  template <typename Body>
  static void run(Body&& body) {
    SEXP continuation = NULL;
    try {
      callConvertingExceptions(Name.text, std::forward<Body>(body));
    } catch (const UnwindJump& jump) {
      continuation = jump.continuation;
    }
    if (continuation != NULL) R_ContinueUnwind(continuation); // does not return
  }
};

} // namespace bartcore_bridge

#endif // R_INTERFACE_EXCEPTIONS_HPP
