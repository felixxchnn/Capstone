// Headless-Chrome smoke test for the built UI.
// Serves ./dist on a local port and drives a real browser through the routes,
// capturing screenshots and console errors. Not a unit test; a sanity check.
//
//   node scripts/browser-smoke.mjs            # desktop 1440x900
//   node scripts/browser-smoke.mjs --mobile   # 390x844
//
// Requires a local Chrome/Edge. Writes screenshots to ./.smoke/ (git-ignored).

import { createServer } from 'node:http';
import { readFile, mkdir, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { extname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { spawn } from 'node:child_process';

const here = fileURLToPath(new URL('.', import.meta.url));
const dist = resolve(here, '..', 'dist');
const outDir = resolve(here, '..', '.smoke');
const mobile = process.argv.includes('--mobile');
const PORT = 4321;

const MIME = {
  '.html': 'text/html',
  '.js': 'text/javascript',
  '.css': 'text/css',
  '.json': 'application/json',
  '.svg': 'image/svg+xml',
  '.jpg': 'image/jpeg',
  '.png': 'image/png',
  '.map': 'application/json',
  '.woff2': 'font/woff2',
};

if (!existsSync(dist)) {
  console.error('dist/ not found — run `npm run build` first.');
  process.exit(1);
}

function findBrowser() {
  const cands = [
    'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',
    'C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe',
    'C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe',
    'C:\\Program Files\\Microsoft\\Edge\\Application\\msedge.exe',
  ];
  return cands.find((c) => existsSync(c)) ?? null;
}

const server = createServer(async (req, res) => {
  try {
    const url = new URL(req.url, `http://localhost:${PORT}`);
    let p = join(dist, decodeURIComponent(url.pathname));
    if (!existsSync(p) || (await stat(p)).isDirectory()) p = join(dist, 'index.html');
    const body = await readFile(p);
    res.writeHead(200, { 'content-type': MIME[extname(p)] ?? 'application/octet-stream' });
    res.end(body);
  } catch (err) {
    res.writeHead(500);
    res.end(String(err));
  }
});

await new Promise((r) => server.listen(PORT, r));
await mkdir(outDir, { recursive: true });

const browser = findBrowser();
if (!browser) {
  console.error('No local Chrome/Edge found; skipping browser smoke.');
  server.close();
  process.exit(0);
}

const routes = [
  ['overview', '/'],
  ['explore', '/#/explore'],
  ['compare', '/#/compare'],
  ['methods', '/#/methods'],
  ['structure', '/#/structure'],
];

const width = mobile ? 390 : 1440;
const height = mobile ? 844 : 900;
let failures = 0;

for (const [name, route] of routes) {
  const shot = join(outDir, `${mobile ? 'mobile' : 'desktop'}-${name}.png`);
  const userDir = join(outDir, `prof-${name}`);
  const args = [
    '--headless=new',
    '--disable-gpu',
    '--no-sandbox',
    '--hide-scrollbars',
    `--user-data-dir=${userDir}`,
    `--window-size=${width},${height}`,
    '--virtual-time-budget=8000',
    '--run-all-compositor-stages-before-draw',
    `--screenshot=${shot}`,
    `http://localhost:${PORT}${route}`,
  ];
  const { code, stderr } = await run(browser, args);
  // Only page-level JS problems count. Chrome writes many unrelated
  // ERROR: lines to stderr in headless mode (USB device log, GCM endpoint,
  // Google-app install, GPU/Vulkan/Fontconfig) — none are our app.
  const consoleErrors = stderr
    .split('\n')
    .filter((l) => /Uncaught|Unhandled|SyntaxError|TypeError:|ReferenceError|\[ErrorBoundary/i.test(l))
    .filter((l) => !/DevTools listening/i.test(l));
  const ok = code === 0 && existsSync(shot) && consoleErrors.length === 0;
  if (!ok) failures += 1;
  console.log(
    `  [${ok ? 'ok' : 'FAIL'}] ${name.padEnd(10)} -> ${shot}` +
      (consoleErrors.length ? `  (${consoleErrors.length} console errors)` : ''),
  );
  if (consoleErrors.length) consoleErrors.slice(0, 4).forEach((l) => console.log('       ' + l.trim()));
}

server.close();
console.log(failures === 0 ? '\nBrowser smoke: pass' : `\nBrowser smoke: ${failures} route(s) failed`);
process.exit(failures === 0 ? 0 : 1);

function run(cmd, args) {
  return new Promise((res) => {
    const p = spawn(cmd, args, { stdio: ['ignore', 'pipe', 'pipe'] });
    let stderr = '';
    p.stderr.on('data', (d) => (stderr += d));
    p.on('close', (code) => res({ code, stderr }));
    setTimeout(() => p.kill(), 30000);
  });
}
