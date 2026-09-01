/** Compact code-native DNA double-helix mark for the nav / footer.
 *  Decorative — aria-hidden; callers provide any needed text. */
export function HelixMark({ size = 28 }: { size?: number }) {
  return (
    <svg
      width={size}
      height={size}
      viewBox="0 0 32 32"
      aria-hidden="true"
      focusable="false"
      className="helix-mark"
    >
      <rect width="32" height="32" rx="7" fill="var(--c-near-black)" />
      <path
        d="M9 6c6 4 8 12 14 20M23 6C17 10 15 18 9 26"
        stroke="var(--c-bio-green)"
        strokeWidth="2.4"
        fill="none"
        strokeLinecap="round"
      />
      <path
        d="M10 10h12M9 16h14M10 22h12"
        stroke="var(--c-mint)"
        strokeWidth="1.6"
        strokeLinecap="round"
      />
    </svg>
  );
}
