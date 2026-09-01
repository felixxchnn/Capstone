/** Detect whether the browser can give us a WebGL2 (or WebGL) context. jsdom
 *  returns null, so tests exercise the text-only fallback path. */
export function detectWebgl(): { ok: boolean; version: 'webgl2' | 'webgl' | null; reason?: string } {
  if (typeof document === 'undefined') {
    return { ok: false, version: null, reason: 'No document (server-side).' };
  }
  let canvas: HTMLCanvasElement;
  try {
    canvas = document.createElement('canvas');
  } catch {
    return { ok: false, version: null, reason: 'Could not create a canvas element.' };
  }
  try {
    if (canvas.getContext('webgl2')) return { ok: true, version: 'webgl2' };
  } catch {
    /* fall through */
  }
  try {
    if (canvas.getContext('webgl') || canvas.getContext('experimental-webgl')) {
      return { ok: true, version: 'webgl', reason: 'Only WebGL 1 is available; rendering may be slower.' };
    }
  } catch {
    /* fall through */
  }
  return {
    ok: false,
    version: null,
    reason:
      'This browser or device reports no WebGL context. Hardware acceleration may be disabled.',
  };
}
