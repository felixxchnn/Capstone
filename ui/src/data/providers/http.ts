// Small fetch helper shared by the structure providers.
// - honours an AbortSignal (stale requests are cancelled when the gene changes)
// - a real timeout
// - distinguishes "offline", "http error", "aborted", "bad json"

export class ProviderError extends Error {
  constructor(
    message: string,
    readonly kind: 'offline' | 'http' | 'timeout' | 'parse' | 'aborted',
    readonly status?: number,
  ) {
    super(message);
    this.name = 'ProviderError';
  }
}

const DEFAULT_TIMEOUT_MS = 12_000;

export async function getJson<T>(
  url: string,
  opts: { signal?: AbortSignal; timeoutMs?: number; accept?: string } = {},
): Promise<T> {
  const res = await rawFetch(url, opts);
  const text = await res.text();
  try {
    return JSON.parse(text) as T;
  } catch {
    throw new ProviderError(`Malformed JSON from ${url}`, 'parse', res.status);
  }
}

export async function getText(
  url: string,
  opts: { signal?: AbortSignal; timeoutMs?: number; accept?: string } = {},
): Promise<string> {
  const res = await rawFetch(url, opts);
  return res.text();
}

async function rawFetch(
  url: string,
  opts: { signal?: AbortSignal; timeoutMs?: number; accept?: string },
): Promise<Response> {
  if (typeof navigator !== 'undefined' && navigator.onLine === false) {
    throw new ProviderError('The browser is offline.', 'offline');
  }
  const timeout = new AbortController();
  const timer = setTimeout(() => timeout.abort(), opts.timeoutMs ?? DEFAULT_TIMEOUT_MS);
  const composite = anySignal([opts.signal, timeout.signal]);
  try {
    const res = await fetch(url, {
      signal: composite,
      headers: opts.accept ? { Accept: opts.accept } : undefined,
      // never send credentials / cookies to structure services
      credentials: 'omit',
      mode: 'cors',
    });
    if (!res.ok) {
      throw new ProviderError(`HTTP ${res.status} from ${url}`, 'http', res.status);
    }
    return res;
  } catch (err) {
    if (err instanceof ProviderError) throw err;
    if (err instanceof DOMException && err.name === 'AbortError') {
      if (opts.signal?.aborted) throw new ProviderError('Request superseded.', 'aborted');
      throw new ProviderError(`Request to ${url} timed out.`, 'timeout');
    }
    throw new ProviderError(
      `Could not reach ${url} (${err instanceof Error ? err.message : 'network error'}).`,
      'offline',
    );
  } finally {
    clearTimeout(timer);
  }
}

/** Combine multiple AbortSignals into one (aborts when any input aborts). */
function anySignal(signals: (AbortSignal | undefined)[]): AbortSignal {
  const controller = new AbortController();
  for (const s of signals) {
    if (!s) continue;
    if (s.aborted) {
      controller.abort();
      break;
    }
    s.addEventListener('abort', () => controller.abort(), { once: true });
  }
  return controller.signal;
}
