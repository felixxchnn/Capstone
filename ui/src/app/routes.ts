export const NAV_ROUTES = [
  {
    path: '/',
    label: 'Overview',
    short: 'Overview',
    description: 'What the project asks and what it found',
  },
  {
    path: '/explore',
    label: 'Dependency Explorer',
    short: 'Explorer',
    description: 'Ranked predicted dependencies, evidence, gene details',
  },
  {
    path: '/compare',
    label: 'Model Comparison',
    short: 'Compare',
    description: 'ridge_pca vs ridge_head — two independent rankings',
  },
  {
    path: '/structure',
    label: 'Protein Structure',
    short: 'Structure',
    description: 'Experimental or predicted structure for the encoded protein',
  },
  {
    path: '/methods',
    label: 'Methods & Limitations',
    short: 'Methods',
    description: 'Pipelines, provenance, hashes, disclaimers',
  },
] as const;

export type NavPath = (typeof NAV_ROUTES)[number]['path'];
