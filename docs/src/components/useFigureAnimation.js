import { onMounted, onBeforeUnmount, watch } from 'vue'

/**
 * Drive a figure's animation with requestAnimationFrame.
 * - `onFrame(dt)` runs once per frame, dt in seconds (capped at 0.05).
 * - The loop runs only while `isPlaying` is true and the figure (`rootRef`) is on screen.
 * - Under `prefers-reduced-motion: reduce`, the figure starts paused; the user can still press play.
 */
export function useFigureAnimation(rootRef, isPlaying, onFrame) {
  let rafId = null
  let lastTs = 0
  let visible = true
  let observer = null

  function frame(ts) {
    rafId = null
    if (!lastTs) lastTs = ts
    const dt = Math.min(0.05, (ts - lastTs) / 1000)
    lastTs = ts
    onFrame(dt)
    schedule()
  }

  function schedule() {
    if (rafId === null && isPlaying.value && visible) rafId = requestAnimationFrame(frame)
  }

  function stop() {
    if (rafId !== null) cancelAnimationFrame(rafId)
    rafId = null
    lastTs = 0
  }

  watch(isPlaying, (playing) => (playing ? schedule() : stop()))

  onMounted(() => {
    if (window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
      isPlaying.value = false
    }
    if ('IntersectionObserver' in window && rootRef.value) {
      observer = new IntersectionObserver((entries) => {
        visible = entries.some((e) => e.isIntersecting)
        if (visible) schedule()
        else stop()
      })
      observer.observe(rootRef.value)
    }
    schedule()
  })

  onBeforeUnmount(() => {
    stop()
    if (observer) observer.disconnect()
  })
}
