/* Copyright (C) 2016-2021 INRAE
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_PARALLEL_FOR_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_PARALLEL_FOR_HPP

#include <atomic>
#include <thread>
#include <vector>

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) ||          \
  defined(_M_IX86)
#include <immintrin.h>
#endif

namespace baryonyx {

inline void
cpu_relax() noexcept
{
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) ||          \
  defined(_M_IX86)
    _mm_pause();
#else
    std::this_thread::yield();
#endif
}

/**
 * Persistent worker threads running a loop body over the indices
 * `[0, count)`. The calling thread takes part in the loop and `run` returns
 * once every index has been processed.
 *
 * The solver runs many loops of a few microseconds each, so a run must cost
 * far less than a thread wake-up: workers busy-wait for a while after each
 * run before sleeping on the generation counter, and the caller busy-waits
 * for the workers to finish.
 */
class parallel_for
{
public:
    /** @param thread_number Total number of threads, the caller included. */
    explicit parallel_for(int thread_number)
    {
        const int workers = thread_number > 1 ? thread_number - 1 : 0;

        m_workers.reserve(static_cast<std::size_t>(workers));
        for (int i = 0; i != workers; ++i)
            m_workers.emplace_back([this] { worker_loop(); });
    }

    parallel_for(const parallel_for&) = delete;
    parallel_for(parallel_for&&) = delete;
    parallel_for& operator=(const parallel_for&) = delete;
    parallel_for& operator=(parallel_for&&) = delete;

    ~parallel_for() noexcept
    {
        m_stop.store(true, std::memory_order_relaxed);
        m_generation.fetch_add(1, std::memory_order_release);
        m_generation.notify_all();
    }

    int thread_number() const noexcept
    {
        return static_cast<int>(m_workers.size()) + 1;
    }

    /** Calls `function(i)` for every `i` in `[0, count)`, each index once. */
    template<typename Function>
    void run(int count, Function& function)
    {
        if (m_workers.empty() || count < 2) {
            for (int i = 0; i < count; ++i)
                function(i);
            return;
        }

        m_body = &function;
        m_invoke = [](void* body, int i) {
            (*static_cast<Function*>(body))(i);
        };
        m_count = count;
        m_next.store(0, std::memory_order_relaxed);
        m_pending.store(static_cast<int>(m_workers.size()),
                        std::memory_order_relaxed);

        m_generation.fetch_add(1, std::memory_order_release);
        m_generation.notify_all();

        work();

        while (m_pending.load(std::memory_order_acquire) != 0)
            cpu_relax();
    }

private:
    static constexpr int spin_before_sleep = 200'000;

    void work() noexcept
    {
        for (int i = m_next.fetch_add(1, std::memory_order_relaxed);
             i < m_count;
             i = m_next.fetch_add(1, std::memory_order_relaxed))
            m_invoke(m_body, i);
    }

    void worker_loop() noexcept
    {
        unsigned seen = 0;

        for (;;) {
            for (int spin = 0;
                 spin != spin_before_sleep &&
                 m_generation.load(std::memory_order_acquire) == seen;
                 ++spin)
                cpu_relax();

            m_generation.wait(seen, std::memory_order_acquire);
            seen = m_generation.load(std::memory_order_acquire);

            if (m_stop.load(std::memory_order_relaxed))
                return;

            work();
            m_pending.fetch_sub(1, std::memory_order_release);
        }
    }

    void* m_body = nullptr;
    void (*m_invoke)(void*, int) = nullptr;
    int m_count = 0;

    std::atomic<int> m_next{ 0 };
    std::atomic<int> m_pending{ 0 };
    std::atomic<unsigned> m_generation{ 0 };
    std::atomic<bool> m_stop{ false };

    // Declared last: the threads are joined before the members they use are
    // destroyed.
    std::vector<std::jthread> m_workers;
};

} // namespace baryonyx

#endif
