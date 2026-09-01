// sync-case-study.mjs
// =====================
// Copy the committed scientific source of truth
//   ../data/processed/case_study.json
// verbatim into the app tree as
//   src/data/case_study.generated.json
// and record its SHA-256 in
//   src/data/case_study.sha256.txt
//
// This is the ONLY way case-study data enters the UI. The generated file is a
// byte-for-byte copy; `src/tests/caseStudyLoad.test.ts` re-verifies the copy's
// SHA-256 against the value pinned in `src/data/case_study.expected-sha256.txt`
// (which mirrors capstone/data-integrity-hashes.md). If the committed JSON ever
// changes, the sync updates the copy and the test tells you the pinned hash
// needs a reviewed bump.
//
// Runs automatically on predev / prebuild / pretest.

import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(here, '..', '..');
const source = resolve(repoRoot, 'data', 'processed', 'case_study.json');
const outDir = resolve(here, '..', 'src', 'data');
const outJson = resolve(outDir, 'case_study.generated.json');
const outSha = resolve(outDir, 'case_study.sha256.txt');

if (!existsSync(source)) {
  console.error(`[sync-case-study] source not found: ${source}`);
  process.exit(1);
}

const raw = readFileSync(source); // Buffer — copied verbatim, no re-serialisation
const sha256 = createHash('sha256').update(raw).digest('hex');

// Validate it is parseable strict JSON and has the expected schema before we
// let it into the app.
let parsed;
try {
  parsed = JSON.parse(raw.toString('utf8'));
} catch (err) {
  console.error(`[sync-case-study] case_study.json is not valid JSON: ${err}`);
  process.exit(1);
}
if (parsed.schema_version !== 'case-study/1') {
  console.error(
    `[sync-case-study] unexpected schema_version ${JSON.stringify(parsed.schema_version)} (want "case-study/1")`,
  );
  process.exit(1);
}

mkdirSync(outDir, { recursive: true });
writeFileSync(outJson, raw); // verbatim bytes
writeFileSync(outSha, sha256 + '\n', 'utf8');

console.log(`[sync-case-study] copied case_study.json (${raw.length} bytes)`);
console.log(`[sync-case-study] sha256 ${sha256}`);
