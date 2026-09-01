/**
 * Code-native DNA double-helix + molecular-node hero illustration. Purely
 * decorative: aria-hidden with a text alternative supplied by the caller.
 * Subtle drift animation is disabled under prefers-reduced-motion (CSS).
 */
export function HeroHelix() {
  const rungs = Array.from({ length: 13 }, (_, i) => i);
  return (
    <svg
      className="hero-helix"
      viewBox="0 0 320 420"
      role="img"
      aria-hidden="true"
      focusable="false"
      preserveAspectRatio="xMidYMid meet"
    >
      <defs>
        <linearGradient id="helix-strand-a" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0" stopColor="#18A866" />
          <stop offset="1" stopColor="#087A45" />
        </linearGradient>
        <linearGradient id="helix-strand-b" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0" stopColor="#DDF7E8" />
          <stop offset="1" stopColor="#18A866" />
        </linearGradient>
      </defs>

      {/* faint grid backdrop */}
      <g className="hero-helix__grid" stroke="#18A866" strokeWidth="0.5" opacity="0.14">
        {Array.from({ length: 9 }, (_, i) => (
          <line key={`v${i}`} x1={i * 40} y1="0" x2={i * 40} y2="420" />
        ))}
        {Array.from({ length: 11 }, (_, i) => (
          <line key={`h${i}`} x1="0" y1={i * 42} x2="320" y2={i * 42} />
        ))}
      </g>

      {/* two sine strands */}
      <path
        className="hero-helix__strand"
        d="M110 10 C 210 60, 10 110, 110 160 C 210 210, 10 260, 110 310 C 210 360, 60 400, 150 410"
        fill="none"
        stroke="url(#helix-strand-a)"
        strokeWidth="4"
        strokeLinecap="round"
      />
      <path
        className="hero-helix__strand hero-helix__strand--b"
        d="M210 10 C 110 60, 310 110, 210 160 C 110 210, 310 260, 210 310 C 110 360, 260 400, 170 410"
        fill="none"
        stroke="url(#helix-strand-b)"
        strokeWidth="4"
        strokeLinecap="round"
      />

      {/* base-pair rungs + molecular nodes */}
      <g className="hero-helix__rungs">
        {rungs.map((i) => {
          const y = 24 + i * 30;
          const phase = Math.sin((i / rungs.length) * Math.PI * 3);
          const x1 = 160 - phase * 66;
          const x2 = 160 + phase * 66;
          return (
            <g key={i}>
              <line x1={x1} y1={y} x2={x2} y2={y} stroke="#087A45" strokeWidth="2" opacity="0.55" />
              <circle cx={x1} cy={y} r="3.4" fill="#18A866" />
              <circle cx={x2} cy={y} r="3.4" fill="#DDF7E8" stroke="#18A866" strokeWidth="1" />
            </g>
          );
        })}
      </g>

      {/* a few free molecular nodes for depth */}
      <g className="hero-helix__free" fill="#18A866">
        <circle cx="40" cy="70" r="4" opacity="0.5" />
        <circle cx="285" cy="150" r="5" opacity="0.4" />
        <circle cx="55" cy="330" r="3" opacity="0.6" />
        <circle cx="270" cy="360" r="4" opacity="0.45" />
      </g>
    </svg>
  );
}
