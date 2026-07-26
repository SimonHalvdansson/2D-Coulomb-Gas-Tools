# Baked shape fields

Each character is stored as a separate `SHP1` binary file named by its
four-digit Unicode code point. For example, `A` is `0041.bin` and `a` is
`0061.bin`.

The 32-byte little-endian header is:

| Offset | Type | Meaning |
|---:|---|---|
| 0 | 4 bytes | `SHP1` magic |
| 4 | `uint16` | Format version |
| 6 | `uint16` | Field width and height |
| 8 | `uint16` | Mask width and height |
| 10 | `uint16` | Quadrature sample count |
| 12 | `float32` | Target area |
| 16 | `float32` | Minimum potential |
| 20 | `float32` | Maximum potential |
| 24 | `uint32` | Number of Float16 field values |
| 28 | `uint32` | Number of packed mask bytes |

The header is followed by interleaved little-endian Float16 values in
`gradientX, gradientY, insideMask, potential` order, then a bit-packed
high-resolution glyph mask.

Regenerate the assets on Windows from the repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File `
  .\tools\render_shape_masks.ps1 -OutputDirectory .\shape_masks_tmp

python .\tools\generate_shape_assets.py `
  --masks .\shape_masks_tmp --output .\shape_assets --workers 4
```
