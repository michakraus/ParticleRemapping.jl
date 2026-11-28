// .vitepress/theme/index.ts
//
// DocumenterVitepress's theme entry without the version picker, the author components and the
// readability menu, which the manual does not use. `custom.css` sets the brand colours and
// content width and is loaded last.
import DefaultTheme from 'vitepress/theme'
import type { Theme as ThemeConfig } from 'vitepress'
import 'virtual:mathjax-styles.css';

// __DV_PLUGIN_COMPONENT_IMPORTS__
import AlgorithmPlayer from '@/AlgorithmPlayer.vue'

import { enhanceAppWithTabs } from 'vitepress-plugin-tabs/client'

import './style.css'
import './docstrings.css' // DocumenterVitepress's default is copied into the build.
import './custom.css'

export const Theme: ThemeConfig = {
  extends: DefaultTheme,
  enhanceApp({ app }) {
    enhanceAppWithTabs(app);
    // __DV_PLUGIN_COMPONENT_REGISTRATIONS__
    app.component('AlgorithmPlayer', AlgorithmPlayer)
  }
}
export default Theme
