<script setup>
import { ref, computed, onMounted, useId } from 'vue'
import { useFigureAnimation } from './useFigureAnimation.js'

const props = defineProps({
  player: {
    type: Object,
    required: true
  }
})

const HOLD = 1.2
const MOVE = 1.6
const STEP_DUR = HOLD + MOVE
const PANEL = { x0: 30, y0: 24, w: 420, h: 300 }

const uid = useId()
const clipId = `ap-clip-${uid}`
const arrowId = `ap-arrow-${uid}`

const root = ref(null)
const mounted = ref(false)
const stepCount = computed(() => props.player.steps.length)
const hasMultipleSteps = computed(() => stepCount.value > 1)
const isPlaying = ref(props.player.steps.length > 1)
const clock = ref(0)

onMounted(() => {
  mounted.value = true
})

useFigureAnimation(root, isPlaying, (dt) => {
  const n = stepCount.value
  if (n <= 1) return
  clock.value = (clock.value + dt) % (n * STEP_DUR)
})

const px = (x) => {
  const [xmin, xmax] = props.player.window
  return PANEL.x0 + ((x - xmin) / (xmax - xmin)) * PANEL.w
}

const py = (v) => {
  const [, , vmin, vmax] = props.player.window
  return PANEL.y0 + ((vmax - v) / (vmax - vmin)) * PANEL.h
}

const frameState = computed(() => {
  const n = stepCount.value
  if (n <= 1) return { fromIdx: 0, toIdx: 0, alpha: 0, activeIdx: 0 }
  const total = n * STEP_DUR
  const t = ((clock.value % total) + total) % total
  const fromIdx = Math.min(n - 1, Math.floor(t / STEP_DUR))
  const local = t - fromIdx * STEP_DUR
  if (local <= HOLD) {
    return { fromIdx, toIdx: fromIdx, alpha: 0, activeIdx: fromIdx }
  }
  const toIdx = (fromIdx + 1) % n
  const alpha = Math.min(1, Math.max(0, (local - HOLD) / MOVE))
  return { fromIdx, toIdx, alpha, activeIdx: toIdx }
})

const activeIndex = computed(() => frameState.value.activeIdx)
const activeStep = computed(() => props.player.steps[activeIndex.value] ?? {
  caption: '',
  x: [],
  v: [],
  r: [],
  c: [],
  targets: null,
  overlays: []
})

const parsedColors = computed(() =>
  props.player.steps.map((st) =>
    st.c.map((hex) => {
      const n = parseInt(hex.slice(1), 16)
      return [(n >> 16) & 255, (n >> 8) & 255, n & 255]
    })
  )
)

const renderedMarkers = computed(() => {
  const { fromIdx, toIdx, alpha } = frameState.value
  const steps = props.player.steps
  if (steps.length === 0) return []
  const s0 = steps[fromIdx]
  const N = props.player.markers
  const out = new Array(N)
  if (alpha === 0 || fromIdx === toIdx) {
    for (let i = 0; i < N; i++) {
      out[i] = {
        cx: px(s0.x[i]).toFixed(2),
        cy: py(s0.v[i]).toFixed(2),
        r: s0.r[i].toFixed(2),
        fill: s0.c[i]
      }
    }
    return out
  }
  const s1 = steps[toIdx]
  const c0 = parsedColors.value[fromIdx]
  const c1 = parsedColors.value[toIdx]
  const w0 = 1 - alpha
  for (let i = 0; i < N; i++) {
    const x = w0 * s0.x[i] + alpha * s1.x[i]
    const v = w0 * s0.v[i] + alpha * s1.v[i]
    const r = w0 * s0.r[i] + alpha * s1.r[i]
    const cr = Math.round(w0 * c0[i][0] + alpha * c1[i][0])
    const cg = Math.round(w0 * c0[i][1] + alpha * c1[i][1])
    const cb = Math.round(w0 * c0[i][2] + alpha * c1[i][2])
    out[i] = {
      cx: px(x).toFixed(2),
      cy: py(v).toFixed(2),
      r: r.toFixed(2),
      fill: `rgb(${cr},${cg},${cb})`
    }
  }
  return out
})

const renderedTargets = computed(() => {
  const tg = activeStep.value.targets
  if (!tg) return []
  const out = new Array(tg.x.length)
  for (let i = 0; i < tg.x.length; i++) {
    out[i] = {
      cx: px(tg.x[i]).toFixed(2),
      cy: py(tg.v[i]).toFixed(2)
    }
  }
  return out
})

const renderedOverlays = computed(() =>
  activeStep.value.overlays.map((o) => {
    if (o.kind === 'line' || o.kind === 'arrow') {
      return {
        kind: o.kind,
        x1: px(o.x1).toFixed(2),
        y1: py(o.v1).toFixed(2),
        x2: px(o.x2).toFixed(2),
        y2: py(o.v2).toFixed(2),
        color: o.color,
        width: o.width,
        dash: o.kind === 'line' && o.dash ? '4,3' : undefined
      }
    }
    if (o.kind === 'rect') {
      const ax = px(o.x)
      const bx = px(o.x + o.dx)
      const ay = py(o.v)
      const by = py(o.v + o.dv)
      return {
        kind: 'rect',
        x: Math.min(ax, bx).toFixed(2),
        y: Math.min(ay, by).toFixed(2),
        width: Math.abs(bx - ax).toFixed(2),
        height: Math.abs(by - ay).toFixed(2),
        color: o.color,
        opacity: o.opacity
      }
    }
    return {
      kind: 'text',
      x: px(o.x).toFixed(2),
      y: py(o.v).toFixed(2),
      text: o.text,
      color: o.color,
      size: o.size
    }
  })
)

function togglePlay() {
  if (!hasMultipleSteps.value) return
  isPlaying.value = !isPlaying.value
}

function stepBy(delta) {
  const n = stepCount.value
  if (n <= 1) return
  isPlaying.value = false
  const nextIdx = (activeIndex.value + delta + n) % n
  clock.value = nextIdx * STEP_DUR
}
</script>

<template>
  <figure
    ref="root"
    class="algorithm-player"
    :data-key="player.id"
    :data-steps="player.steps.length"
  >
    <div v-if="mounted && hasMultipleSteps" class="ap-controls">
      <button
        type="button"
        class="ap-btn"
        :class="{ active: isPlaying }"
        :aria-label="isPlaying ? 'Pause' : 'Play'"
        @click="togglePlay"
      >
        {{ isPlaying ? 'Pause' : 'Play' }}
      </button>
      <div class="ap-btn-group">
        <button
          type="button"
          class="ap-btn"
          aria-label="Previous step"
          @click="stepBy(-1)"
        >
          Previous
        </button>
        <button
          type="button"
          class="ap-btn"
          aria-label="Next step"
          @click="stepBy(1)"
        >
          Next
        </button>
      </div>
    </div>

    <div class="ap-canvas">
      <svg
        xmlns="http://www.w3.org/2000/svg"
        viewBox="0 0 470 354"
        role="img"
        :aria-label="activeStep.caption"
      >
        <defs>
          <marker
            :id="arrowId"
            markerWidth="6"
            markerHeight="6"
            refX="5"
            refY="3"
            orient="auto"
          >
            <path d="M0,0 L6,3 L0,6 z" fill="#808080" />
          </marker>
          <clipPath :id="clipId">
            <rect x="30" y="24" width="420" height="300" />
          </clipPath>
        </defs>

        <rect width="100%" height="100%" fill="white" />
        <rect x="30" y="24" width="420" height="300" fill="white" stroke="#404040" />
        <text
          x="240"
          y="340"
          font-size="12"
          font-style="italic"
          text-anchor="middle"
          fill="#000000"
        >x</text>
        <text
          x="20"
          y="178"
          font-size="12"
          font-style="italic"
          text-anchor="middle"
          fill="#000000"
        >v</text>

        <g :clip-path="`url(#${clipId})`">
          <template v-for="(o, i) in renderedOverlays" :key="i">
            <line
              v-if="o.kind === 'line'"
              :x1="o.x1"
              :y1="o.y1"
              :x2="o.x2"
              :y2="o.y2"
              :stroke="o.color"
              :stroke-width="o.width"
              :stroke-dasharray="o.dash"
            />
            <line
              v-else-if="o.kind === 'arrow'"
              :x1="o.x1"
              :y1="o.y1"
              :x2="o.x2"
              :y2="o.y2"
              :stroke="o.color"
              :stroke-width="o.width"
              :marker-end="`url(#${arrowId})`"
            />
            <rect
              v-else-if="o.kind === 'rect'"
              :x="o.x"
              :y="o.y"
              :width="o.width"
              :height="o.height"
              :fill="o.color"
              :fill-opacity="o.opacity"
            />
            <text
              v-else
              :x="o.x"
              :y="o.y"
              :font-size="o.size"
              :fill="o.color"
              text-anchor="middle"
            >{{ o.text }}</text>
          </template>

          <circle
            v-for="(t, i) in renderedTargets"
            :key="`t-${i}`"
            :cx="t.cx"
            :cy="t.cy"
            r="4.2"
            fill="none"
            stroke="#D55E00"
          />

          <circle
            v-for="(m, i) in renderedMarkers"
            :key="`m-${i}`"
            class="ap-marker"
            :cx="m.cx"
            :cy="m.cy"
            :r="m.r"
            :fill="m.fill"
            fill-opacity="0.85"
          />
        </g>
      </svg>
    </div>

    <figcaption class="ap-caption" :aria-live="isPlaying ? 'off' : 'polite'">
      <span class="ap-step-counter">{{ activeIndex + 1 }} / {{ stepCount }}</span>
      <span class="ap-step-text">{{ activeStep.caption }}</span>
    </figcaption>
  </figure>
</template>
