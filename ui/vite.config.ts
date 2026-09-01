/// <reference types="vitest/config" />
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { fileURLToPath } from 'node:url';

const srcDir = fileURLToPath(new URL('./src', import.meta.url));

// The connected UI reads ONLY the committed case study for its scientific data.
// It never performs model inference. Live protein-structure services (UniProt,
// RCSB, AlphaFold) are contacted at runtime with gene/protein identifiers only.
export default defineConfig({
  plugins: [react()],
  base: './',
  resolve: {
    alias: { '@': srcDir },
  },
  build: {
    target: 'es2022',
    sourcemap: true,
    // Mol* is large; keep it out of the initial bundle. The structures page
    // lazy-imports it, and Rollup will emit it as its own chunk.
    chunkSizeWarningLimit: 4096,
  },
  test: {
    globals: true,
    environment: 'jsdom',
    setupFiles: ['./src/setupTests.ts'],
    include: ['src/tests/**/*.{test,spec}.{ts,tsx}'],
    css: false,
    restoreMocks: true,
    clearMocks: true,
    coverage: {
      provider: 'v8',
      reporter: ['text', 'html'],
      include: ['src/**/*.{ts,tsx}'],
      exclude: ['src/tests/**', 'src/**/*.d.ts', 'src/main.tsx'],
    },
  },
});
