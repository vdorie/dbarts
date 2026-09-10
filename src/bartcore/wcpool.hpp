#ifndef BARTCORE_WCPOOL_HPP
#define BARTCORE_WCPOOL_HPP

#include <atomic>
#include <barrier>
#include <cstddef>
#include <thread>
#include <type_traits>
#include <vector>

/// ThreadSanitizer does not model libc++'s std::barrier as establishing
/// happens-before on this toolchain: a minimal two-thread barrier program
/// reports a race on data the barrier plainly orders. Under the sanitizer the
/// pool therefore states the ordering itself, so a report that survives is
/// about the partitioning rather than about the barrier. The annotations
/// compile to nothing in every other build.
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
#    define BARTCORE_WC_TSAN 1
#  endif
#endif
#ifdef BARTCORE_WC_TSAN
extern "C" void __tsan_acquire(void*);
extern "C" void __tsan_release(void*);
#  define BARTCORE_WC_RELEASE(_P_) __tsan_release(_P_)
#  define BARTCORE_WC_ACQUIRE(_P_) __tsan_acquire(_P_)
#else
#  define BARTCORE_WC_RELEASE(_P_) ((void) 0)
#  define BARTCORE_WC_ACQUIRE(_P_) ((void) 0)
#endif

namespace bartcore {

/// Persistent within-chain worker pool over a std::barrier. The owning thread
/// participates as worker 0, so `participants` helper cores minus one are
/// spawned; workers park on the gate between regions and never call into R.
/// forRange partitions [0, total) into `participants` contiguous slices and
/// runs fn(begin, end) on each, the owner included, blocking until all return.
/// The slice boundaries are a fixed function of total and participants, never
/// of scheduling, which is what lets a reduction built on top stay bitwise
/// invariant across thread count.
class WithinChainPool {
public:
  explicit WithinChainPool(std::size_t participants)
      : n_(participants < 1 ? 1 : participants),
        gate_(static_cast<std::ptrdiff_t>(n_)),
        done_(static_cast<std::ptrdiff_t>(n_)) {
    for (std::size_t w = 1; w < n_; ++w)
      workers_.emplace_back([this, w] { workerLoop(w); });
  }

  ~WithinChainPool() {
    stop_.store(true, std::memory_order_relaxed);
    if (n_ > 1) {
      arriveAt(gate_);
      arriveAt(done_);
    }
    for (std::thread& t : workers_) t.join();
  }

  WithinChainPool(const WithinChainPool&) = delete;
  WithinChainPool& operator=(const WithinChainPool&) = delete;

  std::size_t size() const { return n_; }

  template <class F>
  void forRange(std::size_t total, F&& fn) {
    if (n_ <= 1) {
      fn(std::size_t(0), total);
      return;
    }
    total_ = total;
    JobHolder<std::decay_t<F>> holder{fn};
    job_ = &holder;
    trampoline_ = &JobHolder<std::decay_t<F>>::call;
    arriveAt(gate_);  // open the region
    runRange(0);      // owner runs slice 0
    arriveAt(done_);  // join
    job_ = nullptr;
    trampoline_ = nullptr;
  }

private:
  template <class F>
  struct JobHolder {
    F fn;
    static void call(void* self, std::size_t b, std::size_t e) {
      static_cast<JobHolder*>(self)->fn(b, e);
    }
  };

  static void arriveAt(std::barrier<>& b) {
    BARTCORE_WC_RELEASE(&b);
    b.arrive_and_wait();
    BARTCORE_WC_ACQUIRE(&b);
  }

  void runRange(std::size_t w) {
    std::size_t per = total_ / n_, rem = total_ % n_;
    std::size_t begin = w * per + (w < rem ? w : rem);
    std::size_t count = per + (w < rem ? 1 : 0);
    trampoline_(job_, begin, begin + count);
  }

  void workerLoop(std::size_t w) {
    for (;;) {
      arriveAt(gate_);
      if (stop_.load(std::memory_order_relaxed)) {
        arriveAt(done_);
        return;
      }
      runRange(w);
      arriveAt(done_);
    }
  }

  std::size_t n_;
  std::barrier<> gate_, done_;
  std::vector<std::thread> workers_;
  std::atomic<bool> stop_{false};
  std::size_t total_ = 0;
  void* job_ = nullptr;
  void (*trampoline_)(void*, std::size_t, std::size_t) = nullptr;
};

}  // namespace bartcore

#endif  // BARTCORE_WCPOOL_HPP
