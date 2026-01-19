import * as fs from 'node:fs';
import * as path from 'node:path';
import { Resvg } from '@resvg/resvg-js';
import * as turf from '@turf/turf';
import type {
    Feature,
    FeatureCollection,
    LineString,
    MultiLineString,
    Point,
} from 'geojson';
import proj4 from 'proj4';
import { createElement as h } from 'react';
import satori from 'satori';

// Define projections
const WGS84 = 'EPSG:4326';
const WEB_MERCATOR =
    '+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs';

// California state boundary (simplified) - coordinates in [lon, lat]
const CALIFORNIA_BOUNDARY: number[][] = [
    [-123.233256, 42.006186],
    [-122.378853, 42.011663],
    [-121.037003, 41.995232],
    [-120.001861, 41.995232],
    [-119.996384, 40.264519],
    [-120.001861, 38.999346],
    [-118.71478, 38.101128],
    [-117.498899, 37.21934],
    [-116.540435, 36.501861],
    [-115.85034, 35.970598],
    [-114.634459, 35.00118],
    [-114.634459, 34.87521],
    [-114.470151, 34.710902],
    [-114.333228, 34.448009],
    [-114.136058, 34.305608],
    [-114.256551, 34.174162],
    [-114.415382, 34.108438],
    [-114.535874, 33.933176],
    [-114.497536, 33.697668],
    [-114.524921, 33.54979],
    [-114.727567, 33.40739],
    [-114.661844, 33.034958],
    [-114.524921, 33.029481],
    [-114.470151, 32.843265],
    [-114.524921, 32.755634],
    [-114.72209, 32.717295],
    [-116.04751, 32.624187],
    [-117.126467, 32.536556],
    [-117.24696, 32.668003],
    [-117.252437, 32.876127],
    [-117.329114, 33.122589],
    [-117.471515, 33.297851],
    [-117.7837, 33.538836],
    [-118.183517, 33.763391],
    [-118.260194, 33.703145],
    [-118.413548, 33.741483],
    [-118.391641, 33.840068],
    [-118.566903, 34.042715],
    [-118.802411, 33.998899],
    [-119.218659, 34.146777],
    [-119.278905, 34.26727],
    [-119.558229, 34.415147],
    [-119.875891, 34.40967],
    [-120.138784, 34.475393],
    [-120.472878, 34.448009],
    [-120.64814, 34.579455],
    [-120.609801, 34.858779],
    [-120.670048, 34.902595],
    [-120.631709, 35.099764],
    [-120.894602, 35.247642],
    [-120.905556, 35.450289],
    [-121.004141, 35.461243],
    [-121.168449, 35.636505],
    [-121.283465, 35.674843],
    [-121.332757, 35.784382],
    [-121.716143, 36.195153],
    [-121.896882, 36.315645],
    [-121.935221, 36.638785],
    [-121.858544, 36.6114],
    [-121.787344, 36.803093],
    [-121.929744, 36.978355],
    [-122.105006, 36.956447],
    [-122.335038, 37.115279],
    [-122.417192, 37.241248],
    [-122.400761, 37.361741],
    [-122.515777, 37.520572],
    [-122.515777, 37.783465],
    [-122.329561, 37.783465],
    [-122.406238, 38.15042],
    [-122.488392, 38.112082],
    [-122.504823, 37.931343],
    [-122.701993, 37.893004],
    [-122.937501, 38.029928],
    [-122.97584, 38.265436],
    [-123.129194, 38.451652],
    [-123.331841, 38.566668],
    [-123.44138, 38.698114],
    [-123.737134, 38.95553],
    [-123.687842, 39.032208],
    [-123.824765, 39.366301],
    [-123.764519, 39.552517],
    [-123.85215, 39.831841],
    [-124.109566, 40.10553],
    [-124.361506, 40.259108],
    [-124.410798, 40.439847],
    [-124.158859, 40.877937],
    [-124.109566, 41.025814],
    [-124.158859, 41.14083],
    [-124.065751, 41.442061],
    [-124.147905, 41.715847],
    [-124.257444, 41.781569],
    [-124.213628, 42.000709],
    [-123.233256, 42.006186],
];

// Data source URLs
const DATA_SOURCES = {
    alignment:
        'https://services3.arcgis.com/rGGp0aiv6Rf11t2H/arcgis/rest/services/Alignment_Hybrid_Feb_8_2024/FeatureServer/4/query',
    closures:
        'https://services3.arcgis.com/rGGp0aiv6Rf11t2H/arcgis/rest/services/Closures_and_Construction_Projects_Line_Detail_Read_Only/FeatureServer/3/query',
    completedSections:
        'https://services3.arcgis.com/rGGp0aiv6Rf11t2H/arcgis/rest/services/Closures_and_Detours_Public/FeatureServer/0/query',
    constructionPoints:
        'https://services3.arcgis.com/rGGp0aiv6Rf11t2H/arcgis/rest/services/Closures_and_Detours_Public/FeatureServer/0/query',
};

// Configuration
const CONFIG = {
    alignmentSimplifyTolerance: 0.0015, // keep branch detail in overview map
    // Per-point radius intervals (replaces cluster-based approach)
    // Larger radii needed because points represent project locations, not precise guideway coverage
    completedRadiusKm: 1.2, // Each completed point covers ±1.2 km (increased to capture ~80mi completed guideway)
    constructionRadiusKm: 0.8, // Each construction point covers ±0.8 km
    enableGapFilling: true,
    // Gap-filling: if two snapped points are within this distance, fill the gap
    fillGapKm: 4.0, // Bridge gaps between project markers in same corridor
    gridStepKm: 0.05,
    lineSampleStepKm: 0.1,
    manualCorridorPaddingKm: 0.1,
    manualCorridorSnapKm: 2,
    segmentLengthKm: 1,
    simplifyTolerance: 0.005, // degrees (~500m), for Douglas-Peucker - smoother lines
    snapThresholdKmLine: 1.5, // 1.5km for line features
    // Split snap thresholds: points can be farther from alignment (e.g., on nearby roads)
    snapThresholdKmPoint: 0.75, // 750m for point features
    svgPadding: 20,
    svgStrokeWidth: 4,
    svgWidth: 800,
};

// Manually defined upgraded corridor interval (Caltrain SF to San Jose)
const MANUAL_UPGRADED_CORRIDOR_ENDPOINTS = {
    end: { lat: 37.3297, lon: -121.902, name: 'San Jose Diridon' },
    start: { lat: 37.7765, lon: -122.3942, name: 'San Francisco 4th & King' },
};

// Diagnostics for tracking snapping behavior
interface SnapDiagnostics {
    totalFeatures: number;
    snappedCount: number;
    rejectedCount: number;
    snapDistances: number[]; // distances in km for snapped features
    avgSnapDistanceKm: number;
    maxSnapDistanceKm: number;
}

interface Diagnostics {
    completedPoints: SnapDiagnostics;
    completedLines: SnapDiagnostics & { sampleCount: number };
    constructionPoints: SnapDiagnostics;
    constructionLines: SnapDiagnostics & { sampleCount: number };
    closureLines: SnapDiagnostics & { sampleCount: number };
    intervalsBeforeGapFill: number;
    intervalsAfterGapFill: number;
    gapsFilled: number;
}

interface Meta {
    alignmentFeatures: number;
    closureFeatures: number;
    // Renamed to reflect data-derived nature
    completedCoverageKm: number;
    completedCoverageMiles: number;
    completedSectionFeatures: number;
    constructionCoverageKm: number;
    constructionCoverageMiles: number;
    constructionPointFeatures: number;
    upgradedCorridorKm: number;
    upgradedCorridorMiles: number;
    lastUpdated: string;
    plannedKm: number;
    plannedMiles: number;
    // Updated config fields
    snapThresholdKmPoint: number;
    snapThresholdKmLine: number;
    lineSampleStepKm: number;
    completedRadiusKm: number;
    constructionRadiusKm: number;
    fillGapKm: number;
    enableGapFilling: boolean;
    gridStepKm: number;
    notes: string[];
    sources: string[];
    totalLengthKm: number;
    totalLengthMiles: number;
    diagnostics: Diagnostics;
}

async function fetchGeoJSONWithPagination(
    baseUrl: string,
    whereClause: string = '1=1',
): Promise<FeatureCollection> {
    const allFeatures: Feature[] = [];
    let offset = 0;
    const recordCount = 1000;
    let hasMore = true;

    console.log(`Fetching from ${baseUrl}...`);

    while (hasMore) {
        const params = new URLSearchParams({
            f: 'geojson',
            outFields: '*',
            resultOffset: offset.toString(),
            resultRecordCount: recordCount.toString(),
            returnGeometry: 'true',
            where: whereClause,
        });

        const url = `${baseUrl}?${params.toString()}`;
        const response = await fetch(url);

        if (!response.ok) {
            throw new Error(`Failed to fetch: ${response.statusText}`);
        }

        const data = (await response.json()) as FeatureCollection;

        if (data.features && data.features.length > 0) {
            allFeatures.push(...data.features);
            console.log(
                `  Fetched ${data.features.length} features (total: ${allFeatures.length})`,
            );

            if (data.features.length < recordCount) {
                hasMore = false;
            } else {
                offset += recordCount;
            }
        } else {
            hasMore = false;
        }
    }

    return {
        features: allFeatures,
        type: 'FeatureCollection',
    };
}

function projectToMercator(coords: number[]): number[] {
    return proj4(WGS84, WEB_MERCATOR, coords);
}

function _projectFromMercator(coords: number[]): number[] {
    return proj4(WEB_MERCATOR, WGS84, coords);
}

function distance(a: number[], b: number[]): number {
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    return Math.sqrt(dx * dx + dy * dy);
}

interface MergedLineResult {
    mainLine: Feature<LineString>;
    mercedBranch: Feature<LineString> | null;
    combinedLine: Feature<LineString>; // For snapping purposes (MultiLineString semantics via coords)
    totalLengthKm: number;
    mainLineLengthKm: number;
    mercedBranchLengthKm: number;
}

function mergeLineStrings(
    features: Feature<LineString | MultiLineString>[],
): MergedLineResult {
    // Extract all individual line segments
    const segments: number[][][] = [];

    for (const feature of features) {
        const geom = feature.geometry;
        if (geom.type === 'LineString') {
            segments.push([...geom.coordinates]);
        } else if (geom.type === 'MultiLineString') {
            for (const line of geom.coordinates) {
                segments.push([...line]);
            }
        }
    }

    if (segments.length === 0) {
        const emptyLine = turf.lineString([[0, 0]]);
        return {
            combinedLine: emptyLine,
            mainLine: emptyLine,
            mainLineLengthKm: 0,
            mercedBranch: null,
            mercedBranchLengthKm: 0,
            totalLengthKm: 0,
        };
    }

    // First, dedupe each segment internally to remove back-and-forth within segments
    const cleanedSegments = segments
        .map((seg) => {
            const gridSize = 0.001; // ~100m grid cells
            const seen = new Set<string>();
            const result: number[][] = [];
            for (const coord of seg) {
                const key = `${Math.round(coord[0] / gridSize)},${Math.round(coord[1] / gridSize)}`;
                if (!seen.has(key)) {
                    seen.add(key);
                    result.push(coord);
                }
            }
            return result;
        })
        .filter((seg) => seg.length >= 2);

    // Find the segment closest to SF (northwest terminus of Phase 1)
    // SF 4th & King is at approximately (-122.394, 37.776)
    const SF_COORD = [-122.394, 37.776];
    let startIdx = 0;
    let minDistToSF = Number.POSITIVE_INFINITY;

    for (let i = 0; i < cleanedSegments.length; i++) {
        const seg = cleanedSegments[i];
        // Check both endpoints
        for (const coord of [seg[0], seg[seg.length - 1]]) {
            const dist = distance(coord, SF_COORD);
            if (dist < minDistToSF) {
                minDistToSF = dist;
                startIdx = i;
            }
        }
    }

    // Greedy algorithm to stitch segments together
    // Start from SF to ensure main corridor is included, Merced branch added via gap-free stitching
    const result: number[][] = [...cleanedSegments[startIdx]];
    const used = new Set<number>([startIdx]);

    while (used.size < cleanedSegments.length) {
        const resultEnd = result[result.length - 1];

        let bestIdx = -1;
        let bestDist = Number.POSITIVE_INFINITY;
        let reverseSegment = false;

        // Find the closest unused segment to append (only append, never prepend)
        for (let i = 0; i < cleanedSegments.length; i++) {
            if (used.has(i)) continue;

            const seg = cleanedSegments[i];
            const segStart = seg[0];
            const segEnd = seg[seg.length - 1];

            // Try connecting segment start to result end
            const d1 = distance(resultEnd, segStart);
            if (d1 < bestDist) {
                bestDist = d1;
                bestIdx = i;
                reverseSegment = false;
            }

            // Try connecting segment end to result end (reverse segment)
            const d2 = distance(resultEnd, segEnd);
            if (d2 < bestDist) {
                bestDist = d2;
                bestIdx = i;
                reverseSegment = true;
            }
        }

        if (bestIdx === -1) break;

        used.add(bestIdx);
        let segToAdd = cleanedSegments[bestIdx];
        if (reverseSegment) {
            segToAdd = [...segToAdd].reverse();
        }

        // Skip first point if it's very close to avoid duplicates
        const startIdx =
            distance(result[result.length - 1], segToAdd[0]) < 0.0001 ? 1 : 0;
        result.push(...segToAdd.slice(startIdx));
    }

    // Global deduplication - remove any point that revisits an area
    // Use finer grid to preserve the Merced branch near the junction
    const gridSize = 0.0005; // ~50m grid cells (finer to preserve branch detail)
    const visited = new Map<string, number>(); // grid cell -> index in cleaned array
    const cleaned: number[][] = [];

    for (let i = 0; i < result.length; i++) {
        const coord = result[i];
        const cellKey = `${Math.round(coord[0] / gridSize)},${Math.round(coord[1] / gridSize)}`;

        if (!visited.has(cellKey)) {
            // First time visiting this cell
            visited.set(cellKey, cleaned.length);
            cleaned.push(coord);
        }
        // If we've been here before, skip it entirely (strict dedup)
    }

    // Collect Merced branch coordinates separately before any trimming
    // Merced branch: longitude -120.6 to -120.3, latitude 37.08 to 37.35
    const mercedBranchCoords: number[][] = [];
    for (const seg of cleanedSegments) {
        for (const coord of seg) {
            if (
                coord[0] > -120.6 &&
                coord[0] < -120.3 &&
                coord[1] > 37.08 &&
                coord[1] < 37.35
            ) {
                mercedBranchCoords.push(coord);
            }
        }
    }

    // Dedupe Merced branch
    const mercedDeduped: number[][] = [];
    const mercedVisited = new Set<string>();
    for (const coord of mercedBranchCoords) {
        const key = `${Math.round(coord[0] / gridSize)},${Math.round(coord[1] / gridSize)}`;
        if (!mercedVisited.has(key)) {
            mercedVisited.add(key);
            mercedDeduped.push(coord);
        }
    }
    // Sort Merced branch by latitude (south to north)
    mercedDeduped.sort((a, b) => a[1] - b[1]);

    // Find the southernmost point and trim to avoid backtracking past LA
    let southernmostIdx = 0;
    let lowestLat = cleaned[0][1];

    for (let i = 1; i < cleaned.length; i++) {
        if (cleaned[i][1] < lowestLat) {
            lowestLat = cleaned[i][1];
            southernmostIdx = i;
        }
    }

    // Trim at the southernmost point (LA terminus)
    const mainLine = cleaned.slice(0, southernmostIdx + 1);

    // Check how much of Merced branch is already in mainLine
    const mainLineSet = new Set<string>();
    for (const coord of mainLine) {
        const key = `${Math.round(coord[0] / gridSize)},${Math.round(coord[1] / gridSize)}`;
        mainLineSet.add(key);
    }

    // Add Merced branch points that aren't already in mainLine
    const mercedToAdd: number[][] = [];
    for (const coord of mercedDeduped) {
        const key = `${Math.round(coord[0] / gridSize)},${Math.round(coord[1] / gridSize)}`;
        if (!mainLineSet.has(key)) {
            mercedToAdd.push(coord);
        }
    }

    // Create main line (SF to LA)
    const mainLineFeature = turf.lineString(mainLine);
    const mainLineLengthKm = turf.length(mainLineFeature, {
        units: 'kilometers',
    });

    // Create Merced branch if we have points
    let mercedBranchFeature: Feature<LineString> | null = null;
    let mercedBranchLengthKm = 0;

    if (mercedToAdd.length > 1) {
        mercedBranchFeature = turf.lineString(mercedToAdd);
        mercedBranchLengthKm = turf.length(mercedBranchFeature, {
            units: 'kilometers',
        });
        console.log(
            `  Stitched ${segments.length} segments: main line ${mainLine.length} pts (${mainLineLengthKm.toFixed(1)} km), Merced branch ${mercedToAdd.length} pts (${mercedBranchLengthKm.toFixed(1)} km)`,
        );
    } else {
        console.log(
            `  Stitched ${segments.length} segments into ${mainLine.length} points (${mainLineLengthKm.toFixed(1)} km)`,
        );
    }

    // Combined line for snapping (includes both branches)
    const combinedCoords =
        mercedToAdd.length > 1 ? [...mainLine, ...mercedToAdd] : mainLine;
    const combinedLine = turf.lineString(combinedCoords);

    return {
        combinedLine,
        mainLine: mainLineFeature,
        mainLineLengthKm,
        mercedBranch: mercedBranchFeature,
        mercedBranchLengthKm,
        totalLengthKm: mainLineLengthKm + mercedBranchLengthKm,
    };
}

function extractLineStrings(
    features: Feature<LineString | MultiLineString>[],
): Feature<LineString>[] {
    const lines: Feature<LineString>[] = [];
    for (const feature of features) {
        const geom = feature.geometry;
        if (!geom) continue;
        if (geom.type === 'LineString') {
            lines.push(turf.lineString(geom.coordinates));
        } else if (geom.type === 'MultiLineString') {
            for (const line of geom.coordinates) {
                lines.push(turf.lineString(line));
            }
        }
    }
    return lines;
}

function simplifyLine(
    line: Feature<LineString>,
    tolerance: number,
): Feature<LineString> {
    return turf.simplify(line, {
        highQuality: true,
        tolerance,
    }) as Feature<LineString>;
}

type AlongKm = number;

type Status = 'planned' | 'construction' | 'completed' | 'upgraded_corridor';

type Interval = { startKm: number; endKm: number; status: Status };

type SegmentSlice = {
    segment: Feature<LineString>;
    startKm: number;
    endKm: number;
    midKm: number;
};

function snapPointToAlignmentAlongKm(
    alignment: Feature<LineString>,
    coord: number[],
    maxSnapKm: number,
): { alongKm: AlongKm; snapCoord: number[]; distKm: number } | null {
    const nearest = turf.nearestPointOnLine(alignment, turf.point(coord), {
        units: 'kilometers',
    });
    const distKm =
        typeof nearest.properties?.dist === 'number'
            ? nearest.properties.dist
            : turf.distance(turf.point(coord), nearest, {
                  units: 'kilometers',
              });

    if (distKm > maxSnapKm) {
        return null;
    }

    let alongKm = nearest.properties?.location;
    if (typeof alongKm !== 'number') {
        const start = turf.point(alignment.geometry.coordinates[0]);
        const sliced = turf.lineSlice(start, nearest, alignment);
        alongKm = turf.length(sliced, { units: 'kilometers' });
    }

    return {
        alongKm,
        distKm,
        snapCoord: nearest.geometry.coordinates,
    };
}

function createSegments(
    line: Feature<LineString>,
    segmentLengthKm: number,
): SegmentSlice[] {
    const totalLength = turf.length(line, { units: 'kilometers' });
    const segments: SegmentSlice[] = [];

    let currentDistance = 0;
    while (currentDistance < totalLength) {
        const start = turf.along(line, currentDistance, {
            units: 'kilometers',
        });
        const endDistance = Math.min(
            currentDistance + segmentLengthKm,
            totalLength,
        );
        const end = turf.along(line, endDistance, { units: 'kilometers' });

        // Create a segment by slicing the line
        const segment = turf.lineSlice(start, end, line);
        const midKm = (currentDistance + endDistance) / 2;
        segments.push({
            endKm: endDistance,
            midKm,
            segment: segment as Feature<LineString>,
            startKm: currentDistance,
        });

        currentDistance = endDistance;
    }

    return segments;
}

function clampIntervals(intervals: Interval[], totalKm: number): Interval[] {
    const clamped: Interval[] = [];
    for (const interval of intervals) {
        const startKm = Math.max(0, interval.startKm);
        const endKm = Math.min(totalKm, interval.endKm);
        if (endKm > startKm) {
            clamped.push({ endKm, startKm, status: interval.status });
        }
    }
    return clamped;
}

function paintIntervals(
    intervals: Interval[],
    totalKm: number,
    stepKm: number,
): { bins: Status[]; paintedIntervals: Interval[] } {
    const binCount = Math.max(1, Math.ceil(totalKm / stepKm));
    const bins: Status[] = Array.from({ length: binCount }, () => 'planned');
    const rank: Record<Status, number> = {
        completed: 3,
        construction: 2,
        planned: 0,
        upgraded_corridor: 1,
    };

    for (const interval of intervals) {
        const startIdx = Math.max(0, Math.floor(interval.startKm / stepKm));
        const endIdx = Math.min(binCount, Math.ceil(interval.endKm / stepKm));
        for (let i = startIdx; i < endIdx; i++) {
            if (rank[interval.status] > rank[bins[i]]) {
                bins[i] = interval.status;
            }
        }
    }

    const paintedIntervals: Interval[] = [];
    let currentStatus = bins[0];
    let currentStartIdx = 0;

    for (let i = 1; i < bins.length; i++) {
        if (bins[i] !== currentStatus) {
            paintedIntervals.push({
                endKm: i * stepKm,
                startKm: currentStartIdx * stepKm,
                status: currentStatus,
            });
            currentStatus = bins[i];
            currentStartIdx = i;
        }
    }

    paintedIntervals.push({
        endKm: totalKm,
        startKm: currentStartIdx * stepKm,
        status: currentStatus,
    });

    return { bins, paintedIntervals };
}

function intervalsToTotals(intervals: Interval[]): Record<Status, number> {
    const totals: Record<Status, number> = {
        completed: 0,
        construction: 0,
        planned: 0,
        upgraded_corridor: 0,
    };
    for (const interval of intervals) {
        totals[interval.status] += interval.endKm - interval.startKm;
    }
    return totals;
}

function intervalsFromAlongSamples(
    alongSamples: number[],
    sampleStepKm: number,
    status: Status,
): Interval[] {
    if (alongSamples.length === 0) return [];
    const sorted = [...alongSamples].sort((a, b) => a - b);
    const gapThreshold = sampleStepKm * 1.5;
    const paddingKm = sampleStepKm / 2;
    const intervals: Interval[] = [];

    let start = sorted[0];
    let prev = sorted[0];

    for (let i = 1; i < sorted.length; i++) {
        const current = sorted[i];
        if (current - prev > gapThreshold) {
            intervals.push({
                endKm: prev + paddingKm,
                startKm: start - paddingKm,
                status,
            });
            start = current;
        }
        prev = current;
    }

    intervals.push({
        endKm: prev + paddingKm,
        startKm: start - paddingKm,
        status,
    });

    return intervals;
}

function _sampleLineToAlongKms(
    line: Feature<LineString>,
    alignment: Feature<LineString>,
    sampleStepKm: number,
    maxSnapKm: number,
): number[] {
    const lineLengthKm = turf.length(line, { units: 'kilometers' });
    const alongSamples: number[] = [];

    if (lineLengthKm === 0) return alongSamples;

    for (let d = 0; d <= lineLengthKm; d += sampleStepKm) {
        const point = turf.along(line, d, { units: 'kilometers' });
        const snapped = snapPointToAlignmentAlongKm(
            alignment,
            point.geometry.coordinates,
            maxSnapKm,
        );
        if (snapped) {
            alongSamples.push(snapped.alongKm);
        }
    }

    const endpoints = [
        line.geometry.coordinates[0],
        line.geometry.coordinates[line.geometry.coordinates.length - 1],
    ];
    for (const coord of endpoints) {
        const snapped = snapPointToAlignmentAlongKm(
            alignment,
            coord,
            maxSnapKm,
        );
        if (snapped) {
            alongSamples.push(snapped.alongKm);
        }
    }

    return alongSamples;
}

function buildLineCoverageIntervals(
    lines: Feature<LineString>[],
    alignment: Feature<LineString>,
    sampleStepKm: number,
    maxSnapKm: number,
    status: Status,
): {
    intervals: Interval[];
    diagnostics: SnapDiagnostics & { sampleCount: number };
} {
    const intervals: Interval[] = [];
    const snapDistances: number[] = [];
    let totalSamples = 0;
    let snappedSamples = 0;
    let rejectedSamples = 0;

    for (const line of lines) {
        const lineLengthKm = turf.length(line, { units: 'kilometers' });
        if (lineLengthKm === 0) continue;

        const alongSamples: number[] = [];

        // Sample along the line
        for (let d = 0; d <= lineLengthKm; d += sampleStepKm) {
            totalSamples++;
            const point = turf.along(line, d, { units: 'kilometers' });
            const snapped = snapPointToAlignmentAlongKm(
                alignment,
                point.geometry.coordinates,
                maxSnapKm,
            );
            if (snapped) {
                snappedSamples++;
                snapDistances.push(snapped.distKm);
                alongSamples.push(snapped.alongKm);
            } else {
                rejectedSamples++;
            }
        }

        // Also sample endpoints
        const endpoints = [
            line.geometry.coordinates[0],
            line.geometry.coordinates[line.geometry.coordinates.length - 1],
        ];
        for (const coord of endpoints) {
            totalSamples++;
            const snapped = snapPointToAlignmentAlongKm(
                alignment,
                coord,
                maxSnapKm,
            );
            if (snapped) {
                snappedSamples++;
                snapDistances.push(snapped.distKm);
                alongSamples.push(snapped.alongKm);
            } else {
                rejectedSamples++;
            }
        }

        intervals.push(
            ...intervalsFromAlongSamples(alongSamples, sampleStepKm, status),
        );
    }

    const diagnostics: SnapDiagnostics & { sampleCount: number } = {
        avgSnapDistanceKm:
            snapDistances.length > 0
                ? snapDistances.reduce((a, b) => a + b, 0) /
                  snapDistances.length
                : 0,
        maxSnapDistanceKm:
            snapDistances.length > 0 ? Math.max(...snapDistances) : 0,
        rejectedCount: rejectedSamples,
        sampleCount: totalSamples,
        snapDistances,
        snappedCount: snappedSamples,
        totalFeatures: lines.length,
    };

    return { diagnostics, intervals };
}

/**
 * Build intervals from points using per-point radius approach.
 * Each snapped point creates an interval [alongKm - radiusKm, alongKm + radiusKm].
 * This makes coverage proportional to the number of points, not cluster span.
 */
function buildPointRadiusIntervals(
    points: FeatureCollection<Point>,
    alignment: Feature<LineString>,
    status: Status,
    maxSnapKm: number,
    radiusKm: number,
): { intervals: Interval[]; diagnostics: SnapDiagnostics } {
    const intervals: Interval[] = [];
    const snapDistances: number[] = [];
    let snappedCount = 0;
    let rejectedCount = 0;

    for (const point of points.features) {
        const coords = point.geometry?.coordinates;
        if (!coords) continue;

        const snapped = snapPointToAlignmentAlongKm(
            alignment,
            coords,
            maxSnapKm,
        );

        if (!snapped) {
            rejectedCount++;
            continue;
        }

        snappedCount++;
        snapDistances.push(snapped.distKm);

        // Create interval centered on the snapped point
        intervals.push({
            endKm: snapped.alongKm + radiusKm,
            startKm: snapped.alongKm - radiusKm,
            status,
        });
    }

    const diagnostics: SnapDiagnostics = {
        avgSnapDistanceKm:
            snapDistances.length > 0
                ? snapDistances.reduce((a, b) => a + b, 0) /
                  snapDistances.length
                : 0,
        maxSnapDistanceKm:
            snapDistances.length > 0 ? Math.max(...snapDistances) : 0,
        rejectedCount,
        snapDistances,
        snappedCount,
        totalFeatures: points.features.length,
    };

    return { diagnostics, intervals };
}

/**
 * Apply gap-filling: if two intervals of the same status are within fillGapKm,
 * merge them by filling the gap. This handles sparse but logically connected features.
 */
function applyGapFilling(
    intervals: Interval[],
    fillGapKm: number,
): { intervals: Interval[]; gapsFilled: number } {
    if (intervals.length === 0) return { gapsFilled: 0, intervals: [] };

    // Group intervals by status
    const byStatus = new Map<Status, Interval[]>();
    for (const interval of intervals) {
        const list = byStatus.get(interval.status) || [];
        list.push(interval);
        byStatus.set(interval.status, list);
    }

    let totalGapsFilled = 0;
    const result: Interval[] = [];

    for (const [_status, statusIntervals] of byStatus) {
        // Sort by startKm
        const sorted = [...statusIntervals].sort(
            (a, b) => a.startKm - b.startKm,
        );

        if (sorted.length === 0) continue;

        const merged: Interval[] = [];
        let current = { ...sorted[0] };

        for (let i = 1; i < sorted.length; i++) {
            const next = sorted[i];
            // If gap is small enough, merge
            if (next.startKm - current.endKm <= fillGapKm) {
                if (next.startKm > current.endKm) {
                    totalGapsFilled++;
                }
                current.endKm = Math.max(current.endKm, next.endKm);
            } else {
                merged.push(current);
                current = { ...next };
            }
        }
        merged.push(current);
        result.push(...merged);
    }

    return { gapsFilled: totalGapsFilled, intervals: result };
}

function _statusAtKm(bins: Status[], stepKm: number, km: number): Status {
    if (bins.length === 0) return 'planned';
    const idx = Math.min(bins.length - 1, Math.max(0, Math.floor(km / stepKm)));
    return bins[idx];
}

function statusForRange(
    bins: Status[],
    stepKm: number,
    startKm: number,
    endKm: number,
): Status {
    if (bins.length === 0) return 'planned';
    const rank: Record<Status, number> = {
        completed: 3,
        construction: 2,
        planned: 0,
        upgraded_corridor: 1,
    };
    const startIdx = Math.max(0, Math.floor(startKm / stepKm));
    const endIdx = Math.min(bins.length, Math.ceil(endKm / stepKm));
    let best: Status = 'planned';
    for (let i = startIdx; i < endIdx; i++) {
        if (rank[bins[i]] > rank[best]) {
            best = bins[i];
        }
    }
    return best;
}

function coordsToSVGPath(
    coords: number[][],
    transform: (coord: number[]) => { x: number; y: number },
): string {
    if (coords.length === 0) return '';

    const points = coords.map(transform);
    let d = `M ${points[0].x.toFixed(2)} ${points[0].y.toFixed(2)}`;

    for (let i = 1; i < points.length; i++) {
        d += ` L ${points[i].x.toFixed(2)} ${points[i].y.toFixed(2)}`;
    }

    return d;
}

function generateSVG(
    segments: { coords: number[][]; status: Status }[],
    alignmentLines: Feature<LineString>[],
    width: number,
    padding: number,
    strokeWidth: number,
): string {
    // Include California boundary in bounds calculation for proper framing
    const allCoords = [
        ...segments.flatMap((s) => s.coords),
        ...alignmentLines.flatMap((line) => line.geometry.coordinates),
        ...CALIFORNIA_BOUNDARY,
    ];
    const projectedAll = allCoords.map(projectToMercator);

    // Calculate bounds
    let minX = Number.POSITIVE_INFINITY;
    let maxX = Number.NEGATIVE_INFINITY;
    let minY = Number.POSITIVE_INFINITY;
    let maxY = Number.NEGATIVE_INFINITY;

    for (const [x, y] of projectedAll) {
        minX = Math.min(minX, x);
        maxX = Math.max(maxX, x);
        minY = Math.min(minY, y);
        maxY = Math.max(maxY, y);
    }

    // Calculate scale and dimensions
    const dataWidth = maxX - minX;
    const dataHeight = maxY - minY;
    const availableWidth = width - 2 * padding;
    const scale = availableWidth / dataWidth;
    const height = dataHeight * scale + 2 * padding;

    // Transform function
    const transform = (coord: number[]): { x: number; y: number } => {
        const [lon, lat] = coord;
        const [x, y] = projectToMercator([lon, lat]);
        return {
            x: (x - minX) * scale + padding,
            y: height - ((y - minY) * scale + padding), // Flip Y axis
        };
    };

    // Generate California state background path
    const californiaBgPath = coordsToSVGPath(CALIFORNIA_BOUNDARY, transform);

    // Color mapping
    const colors = {
        completed: '#22c55e', // green
        construction: '#eab308', // yellow
        planned: '#dc2626', // red
        upgraded_corridor: '#3b82f6', // blue
    };

    const classNames = {
        completed: 'completed-section',
        construction: 'active-construction',
        planned: 'base-route',
        upgraded_corridor: 'upgraded-corridor',
    };

    // Merge consecutive segments with the same status, but only if they're geographically contiguous
    const mergedSegments: { coords: number[][]; status: Status }[] = [];
    // Threshold for considering segments contiguous (in degrees, ~1km)
    const CONTIGUITY_THRESHOLD = 0.01;

    for (const segment of segments) {
        if (mergedSegments.length === 0) {
            mergedSegments.push({
                coords: [...segment.coords],
                status: segment.status,
            });
        } else {
            const last = mergedSegments[mergedSegments.length - 1];
            const lastCoord = last.coords[last.coords.length - 1];
            const firstCoord = segment.coords[0];
            const dist = Math.sqrt(
                (lastCoord[0] - firstCoord[0]) ** 2 +
                    (lastCoord[1] - firstCoord[1]) ** 2,
            );

            // Only merge if same status AND geographically contiguous
            if (last.status === segment.status && dist < CONTIGUITY_THRESHOLD) {
                // Same status and contiguous, merge by appending coordinates
                if (dist < 0.0001) {
                    // Very close, skip duplicate first point
                    last.coords.push(...segment.coords.slice(1));
                } else {
                    last.coords.push(...segment.coords);
                }
            } else {
                // Different status OR not contiguous, start new segment
                mergedSegments.push({
                    coords: [...segment.coords],
                    status: segment.status,
                });
            }
        }
    }

    console.log(
        `  Merged ${segments.length} segments into ${mergedSegments.length} paths`,
    );

    // Generate paths for each merged segment
    const paths = mergedSegments
        .map((segment) => {
            const d = coordsToSVGPath(segment.coords, transform);
            if (!d) return '';
            return `<path
    d="${d}"
    fill="none"
    stroke="${colors[segment.status]}"
    stroke-width="${strokeWidth}"
    stroke-linecap="round"
    stroke-linejoin="round"
    class="${classNames[segment.status]}"
  />`;
        })
        .filter((p) => p.length > 0);

    // Build SVG with California background
    const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height.toFixed(0)}" preserveAspectRatio="xMidYMid meet" class="route-map">
  <!-- California state background -->
  <path
    d="${californiaBgPath} Z"
    fill="#18181b"
    stroke="#27272a"
    stroke-width="1"
    class="california-bg"
  />
  ${paths.join('\n  ')}
</svg>`;

    return svg;
}

async function generateOGImage(
    meta: Meta,
    fontData: ArrayBuffer,
    fontBoldData: ArrayBuffer,
): Promise<Buffer> {
    const width = 1200;
    const height = 630;

    // Calculate percentages
    const percentComplete = Math.round(
        (meta.completedCoverageMiles / meta.totalLengthMiles) * 100,
    );
    const percentConstruction = Math.round(
        (meta.constructionCoverageMiles / meta.totalLengthMiles) * 100,
    );
    const percentUpgraded = Math.round(
        (meta.upgradedCorridorMiles / meta.totalLengthMiles) * 100,
    );
    const percentPlanned = Math.round(
        (meta.plannedMiles / meta.totalLengthMiles) * 100,
    );

    // Colors matching Tailwind
    const colors = {
        amber: '#fcd34d', // amber-300
        amberBg: '#92400e', // amber-800
        barBg: '#27272a', // zinc-800
        bg: '#09090b', // zinc-950
        completed: '#22c55e', // green-500
        construction: '#eab308', // yellow-500
        planned: '#dc2626', // red-600
        text: '#ffffff',
        textMuted: '#a1a1aa', // zinc-400
        upgraded: '#3b82f6', // blue-500
    };

    const legendItems = [
        { color: colors.completed, label: `${percentComplete}% Completed` },
        {
            color: colors.construction,
            label: `${percentConstruction}% Under Construction`,
        },
        {
            color: colors.upgraded,
            label: `${percentUpgraded}% Upgraded Corridor`,
        },
        { color: colors.planned, label: `${percentPlanned}% Planned` },
    ];

    // Create the SVG using satori with h()
    const svg = await satori(
        h(
            'div',
            {
                style: {
                    backgroundColor: colors.bg,
                    display: 'flex',
                    flexDirection: 'column',
                    fontFamily: 'Rubik',
                    height: '100%',
                    justifyContent: 'space-between',
                    padding: '60px 80px',
                    width: '100%',
                },
            },
            h(
                'div',
                {
                    style: {
                        display: 'flex',
                        flexDirection: 'column',
                    },
                },
                h(
                    'div',
                    {
                        style: {
                            color: colors.text,
                            fontSize: 32,
                            fontWeight: 700,
                        },
                    },
                    'Are We Rail Yet?',
                ),
                h(
                    'div',
                    {
                        style: {
                            display: 'flex',
                            transform: 'rotate(-5deg)',
                        },
                    },
                    h(
                        'div',
                        {
                            style: {
                                backgroundColor: colors.amberBg,
                                borderRadius: 8,
                                color: colors.amber,
                                fontSize: 56,
                                fontWeight: 700,
                                padding: '12px 32px',
                            },
                        },
                        'Not yet.',
                    ),
                ),
            ),
            h(
                'div',
                {
                    style: {
                        display: 'flex',
                        flexDirection: 'column',
                        gap: 24,
                    },
                },
                h(
                    'div',
                    {
                        style: {
                            backgroundColor: colors.barBg,
                            borderRadius: 24,
                            display: 'flex',
                            height: 48,
                            overflow: 'hidden',
                            width: '100%',
                        },
                    },
                    h('div', {
                        style: {
                            backgroundColor: colors.completed,
                            height: '100%',
                            width: `${percentComplete}%`,
                        },
                    }),
                    h('div', {
                        style: {
                            backgroundColor: colors.construction,
                            height: '100%',
                            width: `${percentConstruction}%`,
                        },
                    }),
                    h('div', {
                        style: {
                            backgroundColor: colors.upgraded,
                            height: '100%',
                            width: `${percentUpgraded}%`,
                        },
                    }),
                    h('div', {
                        style: {
                            backgroundColor: colors.planned,
                            height: '100%',
                            width: `${percentPlanned}%`,
                        },
                    }),
                ),
                h(
                    'div',
                    {
                        style: {
                            display: 'flex',
                            gap: 32,
                        },
                    },
                    ...legendItems.map((item) =>
                        h(
                            'div',
                            {
                                key: item.label,
                                style: {
                                    alignItems: 'center',
                                    display: 'flex',
                                    gap: 8,
                                },
                            },
                            h('div', {
                                style: {
                                    backgroundColor: item.color,
                                    borderRadius: 4,
                                    height: 24,
                                    width: 24,
                                },
                            }),
                            h(
                                'span',
                                {
                                    style: {
                                        color: colors.text,
                                        fontSize: 18,
                                    },
                                },
                                item.label,
                            ),
                        ),
                    ),
                ),
                h(
                    'div',
                    {
                        style: {
                            color: colors.textMuted,
                            fontSize: 24,
                        },
                    },
                    `California High-Speed Rail Phase 1: ${meta.totalLengthMiles} miles from San Francisco to Los Angeles`,
                ),
            ),
            h(
                'div',
                {
                    style: {
                        color: colors.textMuted,
                        fontSize: 20,
                    },
                },
                'arewerailyet.com',
            ),
        ),
        {
            fonts: [
                {
                    data: fontData,
                    name: 'Rubik',
                    style: 'normal',
                    weight: 400,
                },
                {
                    data: fontBoldData,
                    name: 'Rubik',
                    style: 'normal',
                    weight: 700,
                },
            ],
            height,
            width,
        },
    );

    // Convert SVG to PNG using resvg
    const resvg = new Resvg(svg, {
        fitTo: { mode: 'width', value: width },
    });
    const pngData = resvg.render();
    return pngData.asPng();
}

async function main(): Promise<void> {
    console.log('Starting California HSR data generation...\n');

    // Fetch all data
    console.log('=== Fetching Alignment Data (Phase 1 only) ===');
    const alignmentData = await fetchGeoJSONWithPagination(
        DATA_SOURCES.alignment,
        "PHASE = 'Phase 1'",
    );

    console.log('\n=== Fetching Active Closures ===');
    const closuresData = await fetchGeoJSONWithPagination(
        DATA_SOURCES.closures,
        "Active = 'Active'",
    );

    console.log('\n=== Fetching Construction Points ===');
    const constructionData = await fetchGeoJSONWithPagination(
        DATA_SOURCES.constructionPoints,
        "status <> 'Completed'",
    );

    console.log('\n=== Fetching Completed Sections ===');
    const completedData = await fetchGeoJSONWithPagination(
        DATA_SOURCES.completedSections,
        "status = 'Completed'",
    );

    console.log('\n=== Processing Data ===');

    // Merge alignment into single line (Phase 1 only)
    console.log('Merging Phase 1 alignment features...');
    const alignmentFeatures = alignmentData.features.filter(
        (f): f is Feature<LineString | MultiLineString> =>
            f.geometry?.type === 'LineString' ||
            f.geometry?.type === 'MultiLineString',
    );
    const alignmentLines = extractLineStrings(alignmentFeatures).map((line) =>
        simplifyLine(line, CONFIG.alignmentSimplifyTolerance),
    );
    const mergedResult = mergeLineStrings(alignmentFeatures);

    // Simplify lines separately to avoid creating segments across the LA-to-Merced jump
    console.log('Simplifying lines...');
    const simplifiedMainLine = simplifyLine(
        mergedResult.mainLine,
        CONFIG.simplifyTolerance,
    );
    const simplifiedMercedBranch = mergedResult.mercedBranch
        ? simplifyLine(mergedResult.mercedBranch, CONFIG.simplifyTolerance)
        : null;

    // Create a combined simplified line for snapping (coordinates only, no false jump distance)
    const simplifiedLine = simplifiedMercedBranch
        ? turf.lineString([
              ...simplifiedMainLine.geometry.coordinates,
              ...simplifiedMercedBranch.geometry.coordinates,
          ])
        : simplifiedMainLine;

    // Use the calculated total length (main + Merced branch)
    const totalLengthKm = mergedResult.totalLengthKm;
    console.log(
        `  Total length: ${totalLengthKm.toFixed(2)} km (main + Merced branch)`,
    );
    console.log(
        `  Main line points: ${simplifiedMainLine.geometry.coordinates.length}`,
    );
    if (simplifiedMercedBranch) {
        console.log(
            `  Merced branch points: ${simplifiedMercedBranch.geometry.coordinates.length}`,
        );
    }

    // Create segments separately for each branch (no jump segments)
    console.log('Creating segments...');
    const mainSegments = createSegments(
        simplifiedMainLine,
        CONFIG.segmentLengthKm,
    );

    // For Merced branch, we need to offset the startKm/endKm by the main line length
    // so the interval painting works correctly
    const segments = [...mainSegments];
    if (simplifiedMercedBranch) {
        const mercedSegments = createSegments(
            simplifiedMercedBranch,
            CONFIG.segmentLengthKm,
        );
        // Offset Merced branch segments by main line length
        const mainLineLength = mergedResult.mainLineLengthKm;
        for (const seg of mercedSegments) {
            segments.push({
                ...seg,
                endKm: seg.endKm + mainLineLength,
                midKm: seg.midKm + mainLineLength,
                startKm: seg.startKm + mainLineLength,
            });
        }
    }
    console.log(
        `  Created ${segments.length} segments (main: ${mainSegments.length}, Merced: ${segments.length - mainSegments.length})`,
    );

    // Filter closures and points to valid geometries
    const closures: FeatureCollection<LineString | MultiLineString> = {
        features: closuresData.features.filter(
            (f): f is Feature<LineString | MultiLineString> =>
                f.geometry?.type === 'LineString' ||
                f.geometry?.type === 'MultiLineString',
        ),
        type: 'FeatureCollection',
    };

    const constructionPoints: FeatureCollection<Point> = {
        features: constructionData.features.filter(
            (f): f is Feature<Point> => f.geometry?.type === 'Point',
        ),
        type: 'FeatureCollection',
    };

    const completedPoints: FeatureCollection<Point> = {
        features: completedData.features.filter(
            (f): f is Feature<Point> => f.geometry?.type === 'Point',
        ),
        type: 'FeatureCollection',
    };

    const constructionLines = extractLineStrings(
        constructionData.features.filter(
            (f): f is Feature<LineString | MultiLineString> =>
                f.geometry?.type === 'LineString' ||
                f.geometry?.type === 'MultiLineString',
        ),
    );
    const completedLines = extractLineStrings(
        completedData.features.filter(
            (f): f is Feature<LineString | MultiLineString> =>
                f.geometry?.type === 'LineString' ||
                f.geometry?.type === 'MultiLineString',
        ),
    );

    // Build line coverage intervals with diagnostics (using line snap threshold)
    console.log('\n=== Building Coverage Intervals ===');

    const constructionLineResult = buildLineCoverageIntervals(
        constructionLines,
        simplifiedLine,
        CONFIG.lineSampleStepKm,
        CONFIG.snapThresholdKmLine,
        'construction',
    );
    console.log(
        `  Construction lines: ${constructionLineResult.diagnostics.totalFeatures} features, ${constructionLineResult.diagnostics.snappedCount}/${constructionLineResult.diagnostics.sampleCount} samples snapped`,
    );

    const completedLineResult = buildLineCoverageIntervals(
        completedLines,
        simplifiedLine,
        CONFIG.lineSampleStepKm,
        CONFIG.snapThresholdKmLine,
        'completed',
    );
    console.log(
        `  Completed lines: ${completedLineResult.diagnostics.totalFeatures} features, ${completedLineResult.diagnostics.snappedCount}/${completedLineResult.diagnostics.sampleCount} samples snapped`,
    );

    const closureResult = buildLineCoverageIntervals(
        extractLineStrings(closures.features),
        simplifiedLine,
        CONFIG.lineSampleStepKm,
        CONFIG.snapThresholdKmLine,
        'construction',
    );
    console.log(
        `  Closure lines: ${closureResult.diagnostics.totalFeatures} features, ${closureResult.diagnostics.snappedCount}/${closureResult.diagnostics.sampleCount} samples snapped`,
    );

    // Build point radius intervals with diagnostics (using point snap threshold)
    const constructionPointResult = buildPointRadiusIntervals(
        constructionPoints,
        simplifiedLine,
        'construction',
        CONFIG.snapThresholdKmPoint,
        CONFIG.constructionRadiusKm,
    );
    console.log(
        `  Construction points: ${constructionPointResult.diagnostics.totalFeatures} features, ${constructionPointResult.diagnostics.snappedCount} snapped, ${constructionPointResult.diagnostics.rejectedCount} rejected`,
    );
    if (constructionPointResult.diagnostics.snappedCount > 0) {
        console.log(
            `    Avg snap distance: ${constructionPointResult.diagnostics.avgSnapDistanceKm.toFixed(3)} km, max: ${constructionPointResult.diagnostics.maxSnapDistanceKm.toFixed(3)} km`,
        );
    }

    const completedPointResult = buildPointRadiusIntervals(
        completedPoints,
        simplifiedLine,
        'completed',
        CONFIG.snapThresholdKmPoint,
        CONFIG.completedRadiusKm,
    );
    console.log(
        `  Completed points: ${completedPointResult.diagnostics.totalFeatures} features, ${completedPointResult.diagnostics.snappedCount} snapped, ${completedPointResult.diagnostics.rejectedCount} rejected`,
    );
    if (completedPointResult.diagnostics.snappedCount > 0) {
        console.log(
            `    Avg snap distance: ${completedPointResult.diagnostics.avgSnapDistanceKm.toFixed(3)} km, max: ${completedPointResult.diagnostics.maxSnapDistanceKm.toFixed(3)} km`,
        );
    }

    const manualIntervals: Interval[] = [];
    const manualStart = snapPointToAlignmentAlongKm(
        simplifiedLine,
        [
            MANUAL_UPGRADED_CORRIDOR_ENDPOINTS.start.lon,
            MANUAL_UPGRADED_CORRIDOR_ENDPOINTS.start.lat,
        ],
        CONFIG.manualCorridorSnapKm,
    );
    const manualEnd = snapPointToAlignmentAlongKm(
        simplifiedLine,
        [
            MANUAL_UPGRADED_CORRIDOR_ENDPOINTS.end.lon,
            MANUAL_UPGRADED_CORRIDOR_ENDPOINTS.end.lat,
        ],
        CONFIG.manualCorridorSnapKm,
    );

    if (manualStart && manualEnd) {
        let startKm =
            Math.min(manualStart.alongKm, manualEnd.alongKm) -
            CONFIG.manualCorridorPaddingKm;
        let endKm =
            Math.max(manualStart.alongKm, manualEnd.alongKm) +
            CONFIG.manualCorridorPaddingKm;
        if (startKm < CONFIG.segmentLengthKm) {
            startKm = 0;
        }
        if (startKm > 0 && manualStart.alongKm < CONFIG.segmentLengthKm) {
            startKm = 0;
        }
        if (startKm > 0 && manualStart.alongKm <= CONFIG.manualCorridorSnapKm) {
            startKm = 0;
        }
        const grid = CONFIG.gridStepKm;
        startKm = Math.floor(startKm / grid) * grid;
        endKm = Math.ceil(endKm / grid) * grid;
        manualIntervals.push({
            endKm,
            startKm,
            status: 'upgraded_corridor',
        });
    } else {
        console.log(
            '  Warning: manual upgraded corridor endpoints did not snap to alignment.',
        );
    }

    // Combine all intervals
    const rawIntervals = [
        ...constructionLineResult.intervals,
        ...completedLineResult.intervals,
        ...closureResult.intervals,
        ...constructionPointResult.intervals,
        ...completedPointResult.intervals,
        ...manualIntervals,
    ];

    console.log(`\n  Raw intervals before gap-fill: ${rawIntervals.length}`);
    const intervalsBeforeGapFill = rawIntervals.length;

    // Apply gap-filling if enabled
    let processedIntervals = rawIntervals;
    let gapsFilled = 0;
    if (CONFIG.enableGapFilling) {
        const gapFillResult = applyGapFilling(rawIntervals, CONFIG.fillGapKm);
        processedIntervals = gapFillResult.intervals;
        gapsFilled = gapFillResult.gapsFilled;
        console.log(
            `  After gap-filling (${CONFIG.fillGapKm} km threshold): ${processedIntervals.length} intervals, ${gapsFilled} gaps filled`,
        );
    }

    const clampedIntervals = clampIntervals(processedIntervals, totalLengthKm);
    console.log(`  Clamped intervals: ${clampedIntervals.length}`);

    const { bins, paintedIntervals } = paintIntervals(
        clampedIntervals,
        totalLengthKm,
        CONFIG.gridStepKm,
    );

    console.log(`  Painted intervals: ${paintedIntervals.length}`);

    console.log('Classifying segments by alignment chainage...');
    const classifiedSegments: { coords: number[][]; status: Status }[] = [];
    let plannedCount = 0;
    let constructionCount = 0;
    let completedCount = 0;
    let upgradedCount = 0;

    for (const segment of segments) {
        const status = statusForRange(
            bins,
            CONFIG.gridStepKm,
            segment.startKm,
            segment.endKm,
        );
        if (status === 'completed') {
            completedCount++;
        } else if (status === 'construction') {
            constructionCount++;
        } else if (status === 'upgraded_corridor') {
            upgradedCount++;
        } else {
            plannedCount++;
        }
        classifiedSegments.push({
            coords: segment.segment.geometry.coordinates,
            status,
        });
    }

    const totals = intervalsToTotals(paintedIntervals);
    const completedKm = totals.completed;
    const constructionKm = totals.construction;
    const upgradedCorridorKm = totals.upgraded_corridor;
    const plannedKm = Math.max(
        0,
        totalLengthKm - completedKm - constructionKm - upgradedCorridorKm,
    );

    // Convert to miles (1 km = 0.621371 miles)
    const KM_TO_MILES = 0.621371;
    const completedMiles = completedKm * KM_TO_MILES;
    const constructionMiles = constructionKm * KM_TO_MILES;
    const upgradedCorridorMiles = upgradedCorridorKm * KM_TO_MILES;
    const plannedMiles = plannedKm * KM_TO_MILES;

    console.log(
        `  Planned: ${plannedCount} segments (${plannedKm.toFixed(1)} km / ${plannedMiles.toFixed(1)} mi)`,
    );
    console.log(
        `  Construction: ${constructionCount} segments (${constructionKm.toFixed(1)} km / ${constructionMiles.toFixed(1)} mi)`,
    );
    console.log(
        `  Completed: ${completedCount} segments (${completedKm.toFixed(1)} km / ${completedMiles.toFixed(1)} mi)`,
    );
    console.log(
        `  Upgraded corridor: ${upgradedCount} segments (${upgradedCorridorKm.toFixed(1)} km / ${upgradedCorridorMiles.toFixed(1)} mi)`,
    );

    // Generate SVG
    console.log('\n=== Generating SVG ===');
    const svg = generateSVG(
        classifiedSegments,
        alignmentLines,
        CONFIG.svgWidth,
        CONFIG.svgPadding,
        CONFIG.svgStrokeWidth,
    );
    console.log(`  SVG size: ${(svg.length / 1024).toFixed(2)} KB`);

    // Build diagnostics object (omit snapDistances arrays to keep meta.json small)
    const diagnostics: Diagnostics = {
        closureLines: {
            ...closureResult.diagnostics,
            snapDistances: [],
        },
        completedLines: {
            ...completedLineResult.diagnostics,
            snapDistances: [],
        },
        completedPoints: {
            ...completedPointResult.diagnostics,
            snapDistances: [], // Omit raw data for smaller output
        },
        constructionLines: {
            ...constructionLineResult.diagnostics,
            snapDistances: [],
        },
        constructionPoints: {
            ...constructionPointResult.diagnostics,
            snapDistances: [],
        },
        gapsFilled,
        intervalsAfterGapFill: processedIntervals.length,
        intervalsBeforeGapFill,
    };

    // Create meta information
    const meta: Meta = {
        alignmentFeatures: alignmentData.features.length,
        closureFeatures: closuresData.features.length,
        // Renamed to reflect data-derived nature
        completedCoverageKm: Math.round(completedKm * 10) / 10,
        completedCoverageMiles: Math.round(completedMiles * 10) / 10,
        completedRadiusKm: CONFIG.completedRadiusKm,
        completedSectionFeatures: completedData.features.length,
        constructionCoverageKm: Math.round(constructionKm * 10) / 10,
        constructionCoverageMiles: Math.round(constructionMiles * 10) / 10,
        constructionPointFeatures: constructionData.features.length,
        constructionRadiusKm: CONFIG.constructionRadiusKm,
        diagnostics,
        enableGapFilling: CONFIG.enableGapFilling,
        fillGapKm: CONFIG.fillGapKm,
        gridStepKm: CONFIG.gridStepKm,
        lastUpdated: new Date().toISOString(),
        lineSampleStepKm: CONFIG.lineSampleStepKm,
        notes: [
            'IMPORTANT: "Coverage" values are DATA-DERIVED estimates, not official completed guideway mileage.',
            'Coverage is computed by snapping project features to the Phase 1 alignment and measuring chainage intervals.',
            'Each snapped point contributes a fixed radius of coverage (completedRadiusKm / constructionRadiusKm).',
            'Gap-filling merges nearby intervals of the same status if within fillGapKm threshold.',
            'The "Completed" dataset contains closure/detour projects, not direct rail guideway geometry.',
            'Caltrain electrification tracked as upgraded corridor, not counted as completed HSR.',
        ],
        plannedKm: Math.round(plannedKm * 10) / 10,
        plannedMiles: Math.round(plannedMiles * 10) / 10,
        snapThresholdKmLine: CONFIG.snapThresholdKmLine,
        // Updated config fields
        snapThresholdKmPoint: CONFIG.snapThresholdKmPoint,
        sources: Object.values(DATA_SOURCES),
        totalLengthKm: Math.round(totalLengthKm * 100) / 100,
        totalLengthMiles: Math.round(totalLengthKm * KM_TO_MILES * 10) / 10,
        upgradedCorridorKm: Math.round(upgradedCorridorKm * 10) / 10,
        upgradedCorridorMiles: Math.round(upgradedCorridorMiles * 10) / 10,
    };

    // Write outputs
    const generatedDir = path.join(process.cwd(), 'src', 'generated');
    fs.mkdirSync(generatedDir, { recursive: true });

    const svgPath = path.join(generatedDir, 'route.svg');
    fs.writeFileSync(svgPath, svg);
    console.log(`\nWrote SVG to ${svgPath}`);

    const metaPath = path.join(generatedDir, 'meta.json');
    fs.writeFileSync(metaPath, JSON.stringify(meta, null, 2));
    console.log(`Wrote meta to ${metaPath}`);

    // Generate OG image (PNG for social media compatibility)
    console.log('\n=== Generating OG Image ===');
    // Download Rubik font from Google Fonts
    const fontUrl =
        'https://fonts.gstatic.com/s/rubik/v31/iJWZBXyIfDnIV5PNhY1KTN7Z-Yh-B4i1UA.ttf';
    const fontBoldUrl =
        'https://fonts.gstatic.com/s/rubik/v31/iJWZBXyIfDnIV5PNhY1KTN7Z-Yh-4I-1UA.ttf';
    console.log('  Downloading Rubik font...');
    const [fontResponse, fontBoldResponse] = await Promise.all([
        fetch(fontUrl),
        fetch(fontBoldUrl),
    ]);
    const fontData = await fontResponse.arrayBuffer();
    const fontBoldData = await fontBoldResponse.arrayBuffer();
    console.log('  Generating image...');
    const ogImagePng = await generateOGImage(meta, fontData, fontBoldData);
    const ogImagePath = path.join(generatedDir, 'og-image.png');
    fs.writeFileSync(ogImagePath, ogImagePng);
    console.log(`Wrote OG image to ${ogImagePath}`);

    console.log('\n=== Generation Complete ===');
    console.log(`Alignment features (Phase 1): ${meta.alignmentFeatures}`);
    console.log(`Closure features: ${meta.closureFeatures}`);
    console.log(
        `Construction point features: ${meta.constructionPointFeatures}`,
    );
    console.log(`Completed section features: ${meta.completedSectionFeatures}`);
    console.log(`Total route length: ${meta.totalLengthKm} km`);
}

main().catch((error) => {
    console.error('Generation failed:', error);
    process.exit(1);
});
