"""Bake letter-shaped Coulomb-gas fields into compact Float16 binary assets."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import struct
from pathlib import Path

import numpy as np


MASK_SIZE = 384
FIELD_SIZE = 160
QUADRATURE_POINTS = 8192
LETTER_EXTENT = 0.84
FIELD_MIN = -1.5
FIELD_MAX = 1.5
BARRIER_STRENGTH = 4.0
LOG_SOFTENING = 1e-4
CHARACTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
HEADER = struct.Struct("<4sHHHHfffII")


def mask_contains_grid(mask: np.ndarray, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    u = (xs / LETTER_EXTENT + 1.0) * 0.5
    v = (1.0 - ys / LETTER_EXTENT) * 0.5
    valid = (u >= 0.0) & (u < 1.0) & (v >= 0.0) & (v < 1.0)
    px = np.clip((u * MASK_SIZE).astype(np.int32), 0, MASK_SIZE - 1)
    py = np.clip((v * MASK_SIZE).astype(np.int32), 0, MASK_SIZE - 1)
    return valid & mask[py, px]


def distance_penalty(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    step = (FIELD_MAX - FIELD_MIN) / (FIELD_SIZE - 1)
    x_axis = FIELD_MIN + np.arange(FIELD_SIZE, dtype=np.float32) * step
    y_axis = FIELD_MAX - np.arange(FIELD_SIZE, dtype=np.float32) * step
    xs, ys = np.meshgrid(x_axis, y_axis)
    dist = np.full((FIELD_SIZE, FIELD_SIZE), 1e9, dtype=np.float32)
    dist[mask_contains_grid(mask, xs, ys)] = 0.0
    diagonal = np.float32(np.sqrt(2.0))

    for y in range(FIELD_SIZE):
        for x in range(FIELD_SIZE):
            value = dist[y, x]
            if x > 0:
                value = min(value, dist[y, x - 1] + 1.0)
            if y > 0:
                value = min(value, dist[y - 1, x] + 1.0)
            if x > 0 and y > 0:
                value = min(value, dist[y - 1, x - 1] + diagonal)
            if x + 1 < FIELD_SIZE and y > 0:
                value = min(value, dist[y - 1, x + 1] + diagonal)
            dist[y, x] = value

    for y in range(FIELD_SIZE - 1, -1, -1):
        for x in range(FIELD_SIZE - 1, -1, -1):
            value = dist[y, x]
            if x + 1 < FIELD_SIZE:
                value = min(value, dist[y, x + 1] + 1.0)
            if y + 1 < FIELD_SIZE:
                value = min(value, dist[y + 1, x] + 1.0)
            if x + 1 < FIELD_SIZE and y + 1 < FIELD_SIZE:
                value = min(value, dist[y + 1, x + 1] + diagonal)
            if x > 0 and y + 1 < FIELD_SIZE:
                value = min(value, dist[y + 1, x - 1] + diagonal)
            dist[y, x] = value

    phi = BARRIER_STRENGTH * np.square(dist * step)
    return dist, phi.astype(np.float32), step


def quadrature_samples(mask: np.ndarray) -> np.ndarray:
    inside_pixels = np.flatnonzero(mask.reshape(-1))
    if not inside_pixels.size:
        raise ValueError("Mask is empty")
    sample_count = min(QUADRATURE_POINTS, inside_pixels.size)
    stride = inside_pixels.size / sample_count
    selected = inside_pixels[
        np.minimum(
            inside_pixels.size - 1,
            np.floor((np.arange(sample_count) + 0.5) * stride).astype(np.int64),
        )
    ]
    px = selected % MASK_SIZE
    py = selected // MASK_SIZE
    samples = np.empty((sample_count, 2), dtype=np.float32)
    samples[:, 0] = ((px + 0.5) / MASK_SIZE * 2.0 - 1.0) * LETTER_EXTENT
    samples[:, 1] = (1.0 - (py + 0.5) / MASK_SIZE * 2.0) * LETTER_EXTENT
    return samples


def build_field(mask: np.ndarray) -> tuple[np.ndarray, int, float]:
    samples = quadrature_samples(mask)
    sample_count = samples.shape[0]
    dist, phi, step = distance_penalty(mask)
    x_axis = FIELD_MIN + np.arange(FIELD_SIZE, dtype=np.float32) * step
    y_axis = FIELD_MAX - np.arange(FIELD_SIZE, dtype=np.float32) * step
    xs, ys = np.meshgrid(x_axis, y_axis)
    points = np.column_stack((xs.reshape(-1), ys.reshape(-1))).astype(np.float32)
    field = np.empty((points.shape[0], 4), dtype=np.float32)
    phi_flat = phi.reshape(-1)
    dist_flat = dist.reshape(-1)

    block_size = 256
    for start in range(0, points.shape[0], block_size):
        end = min(points.shape[0], start + block_size)
        dx = points[start:end, 0, None] - samples[None, :, 0]
        dy = points[start:end, 1, None] - samples[None, :, 1]
        radius_squared = dx * dx + dy * dy + LOG_SOFTENING
        inverse = 1.0 / radius_squared
        field[start:end, 0] = np.sum(dx * inverse, axis=1) * (2.0 / sample_count)
        field[start:end, 1] = np.sum(dy * inverse, axis=1) * (2.0 / sample_count)
        field[start:end, 3] = (
            np.sum(np.log(radius_squared), axis=1) / sample_count
            + phi_flat[start:end]
        )

    for y in range(FIELD_SIZE):
        for x in range(FIELD_SIZE):
            index = y * FIELD_SIZE + x
            if dist[y, x] <= 0.0:
                continue
            left = phi[y, max(0, x - 1)]
            right = phi[y, min(FIELD_SIZE - 1, x + 1)]
            top = phi[max(0, y - 1), x]
            bottom = phi[min(FIELD_SIZE - 1, y + 1), x]
            dx_span = 2.0 * step if 0 < x < FIELD_SIZE - 1 else step
            dy_span = 2.0 * step if 0 < y < FIELD_SIZE - 1 else step
            field[index, 0] += (right - left) / dx_span
            field[index, 1] += (top - bottom) / dy_span

    field[:, 2] = (dist_flat == 0.0).astype(np.float32)
    pixel_area = (2.0 * LETTER_EXTENT / MASK_SIZE) ** 2
    area = float(np.count_nonzero(mask) * pixel_area)
    return field, sample_count, area


def write_asset(mask_path: Path, output_path: Path) -> dict[str, float | int | str]:
    mask = np.fromfile(mask_path, dtype=np.uint8)
    if mask.size != MASK_SIZE * MASK_SIZE:
        raise ValueError(f"{mask_path} contains {mask.size} bytes; expected {MASK_SIZE * MASK_SIZE}")
    mask = mask.reshape(MASK_SIZE, MASK_SIZE).astype(bool)
    field, sample_count, area = build_field(mask)
    field16 = field.astype("<f2")
    decoded = field16.astype(np.float32)
    min_potential = float(np.min(decoded[:, 3]))
    max_potential = float(np.max(decoded[:, 3]))
    packed_mask = np.packbits(mask.reshape(-1), bitorder="little")
    header = HEADER.pack(
        b"SHP1",
        1,
        FIELD_SIZE,
        MASK_SIZE,
        sample_count,
        area,
        min_potential,
        max_potential,
        field16.size,
        packed_mask.size,
    )
    output_path.write_bytes(header + field16.tobytes() + packed_mask.tobytes())
    error = np.abs(field - decoded)
    return {
        "file": output_path.name,
        "bytes": output_path.stat().st_size,
        "maxGradientError": float(np.max(error[:, :2])),
        "maxPotentialError": float(np.max(error[:, 3])),
    }


def bake_character(character: str, masks: Path, output: Path) -> tuple[str, dict[str, float | int | str]]:
    stem = f"{ord(character):04X}"
    return character, write_asset(
        masks / f"{stem}.mask",
        output / f"{stem}.bin",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masks", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("shape_assets"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    assets = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(bake_character, character, args.masks, args.output): character
            for character in CHARACTERS
        }
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            character, asset = future.result()
            assets[character] = asset
            print(
                f"[{index:02d}/{len(CHARACTERS)}] {character} -> {asset['file']}",
                flush=True,
            )

    manifest = {
        "format": "SHP1",
        "version": 1,
        "encoding": "float16-le",
        "fieldSize": FIELD_SIZE,
        "maskSize": MASK_SIZE,
        "quadraturePoints": QUADRATURE_POINTS,
        "channels": ["gradientX", "gradientY", "insideMask", "potential"],
        "characters": CHARACTERS,
        "totalBytes": sum(asset["bytes"] for asset in assets.values()),
        "assets": {character: assets[character] for character in CHARACTERS},
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(assets)} assets ({manifest['totalBytes']:,} bytes)")


if __name__ == "__main__":
    main()
