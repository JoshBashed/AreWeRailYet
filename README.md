# Are We Rail Yet?

A progress tracker for California High-Speed Rail Phase 1, visualizing completed, under construction, and planned segments.

## Development

```bash
pnpm install
pnpm generate
pnpm dev
```

## How It Works

### The Website

The website is a straightforward Astro static site:

- Static generation at build time
- Tailwind CSS for styling
- A React component for the interactive map (pan/zoom)
- Data imported directly from the generated TypeScript module

### The Generation Script (`scripts/generate.ts`)

> **Warning**: This script is vibe coded. It uses heuristics and hand-tuned parameters to estimate coverage from imperfect public GIS data. The numbers are approximations, not official statistics from CAHSRA.

The script (described by AI):

1. **Fetches GIS data** from CAHSRA's public ArcGIS services:
   - Phase 1 alignment geometry
   - Active closures/construction projects (line features)
   - Construction and completed point markers

2. **Snaps features to the alignment** using `turf.nearestPointOnLine`, converting geographic positions to chainage (distance along the route in km)

3. **Creates coverage intervals** for each snapped feature:
   - Point features: Each point creates an interval of ±0.8 km around its snapped position
   - Line features: Sampled every 0.1 km, converted to intervals

4. **Gap-fills** nearby intervals of the same status (within 4 km) to handle sparse data

5. **Paints bins** along the alignment with precedence: planned < upgraded_corridor < construction < completed

6. **Generates outputs**: SVG map, `meta.json` with statistics, TypeScript module for imports

## License

MIT ([LICENSE.md](LICENSE.md))
